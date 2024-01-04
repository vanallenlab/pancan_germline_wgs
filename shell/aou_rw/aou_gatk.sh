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
export BOTO_CONFIG=~/code/refs/dotfiles/dfci-g2c.aou-rw.boto.cfg
export BOTO_PATH=~/code/refs/dotfiles/

# Prep working directory structure
for dir in data data/cram_paths cromshell cromshell/inputs cromshell/job_ids \
           cromshell/progress; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done

# Copy necessary code to local disk
gsutil -m cp -r $WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}

# Source .bashrc and set boto config
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
       [ $status == "unknown" ]; then
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

# Check progress of each sample and stage completed samples
while read cancer; do
  n=$( cat sample_lists/$cancer.samples.list | wc -l )
  k=0
  while read sid cram crai; do
    ((k++))
    echo -e "Checking $sid; sample $k of $n $cancer patients"
    code/scripts/check_aou_gatk_hc.py \
      --sample-id "$sid" \
      --status-tsv cromshell/progress/$cancer.gatk_hc.sample_progress.tsv \
      --update-status \
      --unsafe
  done < data/cram_paths/$cancer.cram_paths.tsv
done < cancers.list
while read cancer; do
  cut -f2 cromshell/progress/$cancer.gatk_hc.sample_progress.tsv \
  | sort | uniq -c | sort -nrk1,1 \
  | awk -v cancer=$cancer -v OFS="\t" '{ print cancer, "gatk-hc", $2, $1 }'
done < cancers.list

# Clean up garbage
gsutil -m cp \
  uris_to_delete.list \
  $WORKSPACE_BUCKET/dumpster/dfci_g2c.aou_rw.$( date '+%m_%d_%Y.%Hh%Mm%Ss' ).garbage
rm uris_to_delete.list
# TODO: submit cromwell job that deletes these files using a high-thread / low-mem cloud instance

# TODO: reblock & reheader gVCFs


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
