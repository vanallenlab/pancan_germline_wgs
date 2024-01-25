#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-HC and GATK-SV (module 01) on all AoU samples selected for inclusion in G2C pilot

# Note that this code is designed to be run inside the AoU Researcher Workbench
# See aou_bash_utils.sh for custom function definitions used below


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

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

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

# Main loop for the submission and tracking of GATK-HC workflows
# Before launching, this code will check the status of each sample to determine eligibility
# If already successful, it will stage final outputs and clean up temporary files
# It will also update counts in the main sample progress tracker table
while read cancer; do
  check_status gatk-hc $cancer
  update_status_table gatk-hc
  cleanup_garbage
  submit_workflows gatk-hc $cancer
done < cancers.list

# Main loop for the submission and tracking of GATK-HC gVCF postprocessing
# Includes both reblocking and reheadering
# Before launching, this code will check the status of each sample to determine eligibility
# If already successful, it will stage final outputs and clean up temporary files
# It will also update counts in the main sample progress tracker table
while read cancer; do
  check_status gvcf-pp $cancer
  update_status_table gvcf-pp
  cleanup_garbage
  submit_workflows gvcf-pp $cancer
done < cancers.list


###########
# GATK-SV #
###########

# Main loop for the submission and tracking of GATK-SV workflows
# Before launching, this code will check the status of each sample to determine eligibility
# If already successful, it will stage final outputs and clean up temporary files
# It will also update counts in the main sample progress tracker table
while read cancer; do
  check_status gatk-sv $cancer
  update_status_table gatk-sv
  cleanup_garbage
  submit_workflows gatk-sv $cancer
done < cancers.list

