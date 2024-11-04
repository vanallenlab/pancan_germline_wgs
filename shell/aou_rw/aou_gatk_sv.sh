#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-SV cohort mode pipeline on G2C phase 1

# Note that this code is designed to be run inside the AoU Researcher Workbench
# See gatksv_bash_utils.sh for custom function definitions used below


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in cromshell cromshell/inputs cromshell/job_ids cromshell/progress; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}

# Source .bashrc and bash utility functions
. ~/code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh
. code/refs/gatksv_bash_utils.sh

# Infer workspace number and save as environment variable
export WN=$( get_workspace_number )

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Create dependencies .zip for all GATK-SV module submissions
cd code/wdl/gatk-sv && \
zip gatksv.dependencies.zip *.wdl && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Copy sample & batch information
gsutil -m cp -r \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/batch_info \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
  ~/


##################
# 03 | TrainGCNV #
##################

for iter in 1; do
  while read bid; do
    check_batch_module $bid 03
    # TODO: check status of batch + workflow
    # TODO: Update tracker
    # If not started, failed, or doomed launch:
    submit_batch_module $bid 03
    # TODO: If Succeeded but not staged, stage
    # If Running, do nothing
  done < batch_info/dfci-g2c.gatk-sv.batches.w$WN.list
  # TODO: print tracker to screen
  sleep 60m
done
