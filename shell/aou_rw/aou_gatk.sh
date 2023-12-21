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

# Source .bashrc
. code/dotfiles/aou.rw.bashrc

# Copy sample info to local disk
gsutil -m cp -r $WORKSPACE_BUCKET/data/sample_info/sample_lists ./

# Copy manifest of CRAM/CRAI paths
gsutil -u $GPROJECT -m cp \
  gs://fc-aou-datasets-controlled/v7/wgs/cram/manifest.csv \
  ./cram_manifest.csv

# Make .tsv mapping person_id, cram path, and crai path for each cancer type
for cancer in pancreas; do
  fgrep -wf sample_lists/$cancer.samples.list cram_manifest.csv \
  | sed 's/,/\t/g' | sort -Vk1,1 \
  > data/cram_paths/$cancer.cram_paths.tsv
done


###########
# GATK-HC #
###########

# Launch one GATK-HC workflow for each sample
for cancer in pancreas; do
  if ! [ -e cromshell/progress/$cancer.gatk_hc.sample_progress.tsv ]; then
    touch cromshell/progress/$cancer.gatk_hc.sample_progress.tsv
  fi
  while read sid cram crai; do
    # Check if sample has not been launched or has failed/aborted
    status=$( awk -v FS="\t" -v sid=$sid '{ if ($1==sid) print $2 }' \
              cromshell/progress/$cancer.gatk_hc.sample_progress.tsv )
    if [ -z $status ] || [ $status == "failed" ] || [ $status == "aborted" ]; then
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
done

# Check progress of each sample
for cancer in pancreas; do
  while read sid cram crai; do
    echo $sid
    if ! [ -e cromshell/job_ids/$sid.gatk_hc.job_ids.list ]; then
      echo "not_started"
    else
      jid=$( tail -n1 cromshell/job_ids/$sid.gatk_hc.job_ids.list )
      cromshell-alpha status $jid 2>/dev/null \
      | tail -n2 | jq .status | tr -d '"' | tr '[A-Z]' '[a-z]'
    fi
  done < data/cram_paths/$cancer.cram_paths.tsv  | paste - - \
  > cromshell/progress/$cancer.gatk_hc.sample_progress.tsv
done


###########
# GATK-SV #
###########

# Zip all WDLs into dependencies package
cd code/wdl/gatk-sv && \
zip gatksv.dependencies.zip *.wdl && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Launch one GATK-SV workflow for each sample
for cancer in pancreas; do
  if ! [ -e cromshell/progress/$cancer.gatk_sv.sample_progress.tsv ]; then
    touch cromshell/progress/$cancer.gatk_sv.sample_progress.tsv
  fi
  while read sid cram crai; do
    # Check if sample has not been launched or has failed/aborted
    status=$( awk -v FS="\t" -v sid=$sid '{ if ($1==sid) print $2 }' \
              cromshell/progress/$cancer.gatk_sv.sample_progress.tsv )
    if [ -z $status ] || [ $status == "failed" ] || [ $status == "aborted" ]; then
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
done

# Check progress of each sample
for cancer in pancreas; do
  while read sid cram crai; do
    echo $sid
    if ! [ -e cromshell/job_ids/$sid.gatk_sv.job_ids.list ]; then
      echo "not_started"
    else
      jid=$( tail -n1 cromshell/job_ids/$sid.gatk_sv.job_ids.list )
      cromshell-alpha status $jid 2>/dev/null \
      | tail -n2 | jq .status | tr -d '"' | tr '[A-Z]' '[a-z]'
    fi
  done < data/cram_paths/$cancer.cram_paths.tsv  | paste - - \
  > cromshell/progress/$cancer.gatk_sv.sample_progress.tsv
done


