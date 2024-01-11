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
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in data data/cram_paths cromshell cromshell/inputs cromshell/job_ids \
           cromshell/progress; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}

# Source .bashrc and bash utility functions
. ~/code/refs/dotfiles/aou.rw.bashrc
. code/refs/aou_bash_utils.sh

# Copy sample info to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/data/sample_info/sample_lists ./

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
colorectal
kidney
melanoma
prostate
ovary
EOF
# Note: as of 1/8/24, single-sample processing for cancers were divided among 
# multiple workspaces as follows:
# Main workspace: pancreas, esophagus, stomach, lung, liver, kidney, ovary
# Second workspace: colorectal
# Third workspace: melanoma
# Fourth workspace: prostate

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
  while read sid CRAM CRAI; do
    # Check if sample has not been launched or has failed/aborted
    status=$( awk -v FS="\t" -v sid=$sid '{ if ($1==sid) print $2 }' \
              cromshell/progress/$cancer.gatk_hc.sample_progress.tsv )
    if [ -z $status ] || [ $status == "not_started" ] || \
       [ $status == "failed" ] || [ $status == "aborted" ] || \
       [ $status == "unknown" ]; then
      # Format sample-specific input .json
      echo -e "Submitting GATK-HC for $sid ($cancer cancer)"

      eval "cat << EOF
              $(<code/refs/json/aou.gatk_hc.inputs.template.json)
EOF"  | sed 's/\t//g' | paste -s -d\ \
      > cromshell/inputs/$sid.gatkhc.inputs.json
      # Submit job and add job ID to list of jobs for this sample
      cromshell --no_turtle -t 120 -mc submit \
        --options-json code/refs/json/aou.cromwell_options.default.json \
        code/wdl/gatk-hc/haplotypecaller-gvcf-gatk4.wdl \
        cromshell/inputs/$sid.gatkhc.inputs.json \
      | jq .id | tr -d '"' \
      >> cromshell/job_ids/$sid.gatk_hc.job_ids.list
    fi
  done < data/cram_paths/$cancer.cram_paths.tsv
done < cancers.list

# Check progress of each sample and stage completed samples
# (Function defined in function defined in code/refs/aou_bash_utils.sh)
while read cancer; do
  check_status gatk-hc $cancer
done < cancers.list

# Print table of sample progress
# (Function defined in function defined in code/refs/aou_bash_utils.sh)
update_status_table gatk-hc

# Clean up garbage
# (Function defined in code/refs/aou_bash_utils.sh)
cleanup_garbage

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
  while read sid CRAM CRAI; do
    # Check if sample has not been launched or has failed/aborted
    status=$( awk -v FS="\t" -v sid=$sid '{ if ($1==sid) print $2 }' \
              cromshell/progress/$cancer.gatk_sv.sample_progress.tsv )
    if [ -z $status ] || [ $status == "not_started" ] || \
       [ $status == "failed" ] || [ $status == "aborted" ] || \
       [ $status == "unknown" ]; then
      # Format sample-specific input .json
      echo -e "Submitting GATK-SV for $sid ($cancer cancer)"

      eval "cat << EOF
              $(<code/refs/json/aou.gatk_sv_module_01.inputs.template.json)
EOF"  | sed 's/\t//g' | paste -s -d\ \
      > cromshell/inputs/$sid.gatksv.inputs.json
      # Submit job and add job ID to list of jobs for this sample
      cromshell --no_turtle -t 120 -mc submit \
        --options-json code/refs/json/aou.cromwell_options.default.json \
        --dependencies-zip gatksv.dependencies.zip \
        code/wdl/gatk-sv/GatherSampleEvidence.wdl \
        cromshell/inputs/$sid.gatksv.inputs.json \
      | tail -n4 | jq .id | tr -d '"' \
      >> cromshell/job_ids/$sid.gatk_sv.job_ids.list
    fi
  done < data/cram_paths/$cancer.cram_paths.tsv
done < cancers.list

# Check progress of each sample and stage completed samples
# (Function defined in function defined in code/refs/aou_bash_utils.sh)
while read cancer; do
  check_status gatk-sv $cancer
done < cancers.list

# Print table of sample progress
# (Function defined in function defined in code/refs/aou_bash_utils.sh)
update_status_table gatk-sv

# Clean up garbage
# (Function defined in code/refs/aou_bash_utils.sh)
cleanup_garbage

