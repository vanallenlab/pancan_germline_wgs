#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-HC and GATK-SV (module 01) on all AoU samples selected for inclusion in G2C pilot

# Note that this code is designed to be run inside the AoU Researcher Workbench


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"

# Prep working directory structure
for dir in data data/cram_paths cromshell cromshell/inputs cromshell/job_ids \
           cromshell/progress; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done

# Copy necessary code to local disk
gsutil -m cp -r $WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}

# Source .bashrc
. code/refs/dotfiles/aou.rw.bashrc

# Copy sample info to local disk
gsutil -m cp -r $WORKSPACE_BUCKET/data/sample_info/sample_lists ./

# Copy manifest of CRAM/CRAI paths
gsutil -u $GPROJECT -m cp \
  gs://fc-aou-datasets-controlled/v7/wgs/cram/manifest.csv \
  ./cram_manifest.csv

# Write list of cancers ready to be processed
cat << EOF > cancers.list
pancreas
esophagus
stomach
lung
liver
EOF

# Make .tsv mapping person_id, cram path, and crai path for each cancer type
while read cancer; do
  fgrep -wf sample_lists/$cancer.samples.list cram_manifest.csv \
  | sed 's/,/\t/g' | sort -Vk1,1 \
  > data/cram_paths/$cancer.cram_paths.tsv
done < cancers.list


###########
# GATK-HC #
###########

# Launch one GATK-HC workflow for each sample
while read cancer; do
  if ! [ -e cromshell/progress/$cancer.gatk_hc.sample_progress.tsv ]; then
    touch cromshell/progress/$cancer.gatk_hc.sample_progress.tsv
  fi
  while read sid cram crai; do
    # Check if sample has not been launched or has failed/aborted
    status=$( awk -v FS="\t" -v sid=$sid '{ if ($1==sid) print $2 }' \
              cromshell/progress/$cancer.gatk_hc.sample_progress.tsv )
    if [ -z $status ] || [ $status == "not_started" ] || \
       [ $status == "failed" ] || [ $status == "aborted" ] || \
       [ $status == "staged" ]; then
      # Format sample-specific input .json
      CRAM=$cram CRAI=$crai \
        envsubst < code/refs/json/aou.gatk_hc.inputs.template.json \
      | sed 's/\t//g' | paste -s -d\ \
      > cromshell/inputs/$sid.gatkhc.inputs.json
      # Submit job and add job ID to list of jobs for this sample
      cromshell-alpha submit \
        --options-json code/refs/json/aou.cromwell_options.default.json \
        code/wdl/gatk-hc/haplotypecaller-gvcf-gatk4.wdl \
        cromshell/inputs/$sid.gatkhc.inputs.json \
      | tail -n4 | jq .id | tr -d '"' \
      >> cromshell/job_ids/$sid.gatk_hc.job_ids.list
    fi
  done < data/cram_paths/$cancer.cram_paths.tsv
done < cancers.list

# Check progress of each sample
while read cancer; do
  while read sid cram crai; do
    code/scripts/check_aou_gatk_hc.py \
      --sample-id "$sid" \
      --status-tsv cromshell/progress/$cancer.gatk_hc.sample_progress.tsv \
      --update-status
  done < data/cram_paths/$cancer.cram_paths.tsv
  # Print summary of status to terminal
  cut -f2 cromshell/progress/$cancer.gatk_hc.sample_progress.tsv \
  | sort | uniq -c | sort -nrk1,1 | awk -v OFS="\t" '{ print $2, $1 }'
done < cancers.list

# Clean up finished samples
while read cancer; do
  # Reset lists of files to relocate and delete
  for suf in to_move to_delete; do
    if [ -e cromshell/progress/$cancer.gatk_hc.$suf.list ]; then
      rm cromshell/progress/$cancer.gatk_hc.$suf.list
    fi
  done
  while read sid status; do
    # Only process samples that are reported as succeeded by Cromwell
    if [ $status == "succeeded" ]; then
      most_recent=$( tail -n1 cromshell/job_ids/$sid.gatk_hc.job_ids.list )
      echo "$WORKSPACE_BUCKET/cromwell/outputs/HaplotypeCallerGvcf_GATK4/$most_recent/call-MergeGVCFs/**.g.vcf.gz*" \
      >> cromshell/progress/$cancer.gatk_hc.to_move.list
      while read subid; do
        echo "$WORKSPACE_BUCKET/cromwell/execution/HaplotypeCallerGvcf_GATK4/$subid/**" \
        >> cromshell/progress/$cancer.gatk_hc.to_delete.list
        echo "$WORKSPACE_BUCKET/cromwell/outputs/HaplotypeCallerGvcf_GATK4/$subid/**" \
        >> cromshell/progress/$cancer.gatk_hc.to_delete.list
      done < cromshell/job_ids/$sid.gatk_hc.job_ids.list
      status=staged
    fi
    echo -e "$sid\t$status"
  done < cromshell/progress/$cancer.gatk_hc.sample_progress.tsv \
  > cromshell/progress/$cancer.gatk_hc.sample_progress.tsv2
  # Copy final gVCFs and indexes
  cat cromshell/progress/$cancer.gatk_hc.to_move.list \
  | gsutil -m mv -I $WORKSPACE_BUCKET/dfci-g2c-inputs/aou/gatk-hc/
  # Clean up all intermediate files
  cat cromshell/progress/$cancer.gatk_hc.to_delete.list | gsutil -m rm -I
  # Update sample progress manifest
  mv \
    cromshell/progress/$cancer.gatk_hc.sample_progress.tsv2 \
    cromshell/progress/$cancer.gatk_hc.sample_progress.tsv
done < cancers.list

# TODO: Rename gVCF outputs to remove the "wgs_" prefix from each file

# TODO: reblock & rename gVCFs


###########
# GATK-SV #
###########

# Zip all WDLs into dependencies package
cd code/wdl/gatk-sv && \
zip gatksv.dependencies.zip *.wdl && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Launch one GATK-SV workflow for each sample
while read cancer; do
  if ! [ -e cromshell/progress/$cancer.gatk_sv.sample_progress.tsv ]; then
    touch cromshell/progress/$cancer.gatk_sv.sample_progress.tsv
  fi
  while read sid cram crai; do
    # Check if sample has not been launched or has failed/aborted
    status=$( awk -v FS="\t" -v sid=$sid '{ if ($1==sid) print $2 }' \
              cromshell/progress/$cancer.gatk_sv.sample_progress.tsv )
    if [ -z $status ] || [ $status == "not_started" ] || \
       [ $status == "failed" ] || [ $status == "aborted" ] || \
       [ $status == "staged" ]; then
      # Format sample-specific input .json
      SID=$sid CRAM=$cram CRAI=$crai \
        envsubst < code/refs/json/aou.gatk_sv_module_01.inputs.template.json \
      | sed 's/\t//g' | paste -s -d\ \
      > cromshell/inputs/$sid.gatksv.inputs.json
      # Submit job and add job ID to list of jobs for this sample
      cromshell-alpha submit \
        --options-json code/refs/json/aou.cromwell_options.default.json \
        --dependencies-zip gatksv.dependencies.zip \
        code/wdl/gatk-sv/GatherSampleEvidence.wdl \
        cromshell/inputs/$sid.gatksv.inputs.json \
      | tail -n4 | jq .id | tr -d '"' \
      >> cromshell/job_ids/$sid.gatk_sv.job_ids.list
    fi
  done < data/cram_paths/$cancer.cram_paths.tsv
done < cancers.list

# Check progress of each sample
while read cancer; do
  while read sid cram crai; do

    # TODO: update to new 

    echo $sid
    if ! [ -e cromshell/job_ids/$sid.gatk_sv.job_ids.list ]; then
      echo "not_started"
    else
      status=$( awk -v FS="\t" -v sid=$sid '{ if ($1==sid) print $2 }' \
                cromshell/progress/$cancer.gatk_hc.sample_progress.tsv )
      if [ $status == "succeeded" ] || [ $status == "staged" ]; then
        echo $status
      else
        jid=$( tail -n1 cromshell/job_ids/$sid.gatk_sv.job_ids.list )
        cromshell-alpha status $jid 2>/dev/null \
        | tail -n2 | jq .status | tr -d '"' | tr '[A-Z]' '[a-z]'
      fi
    fi
  done < data/cram_paths/$cancer.cram_paths.tsv  | paste - - \
  > cromshell/progress/$cancer.gatk_sv.sample_progress.tsv
  cut -f2 cromshell/progress/$cancer.gatk_sv.sample_progress.tsv \
  | sort | uniq -c | sort -nrk1,1 | awk -v OFS="\t" '{ print $2, $1 }'
done < cancers.list

# # Clean up finished samples
# TODO: implement this
