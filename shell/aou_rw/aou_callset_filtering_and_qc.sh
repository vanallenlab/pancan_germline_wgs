#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Quality control and filtering of G2C germline callset after joint genotyping

# Note that this code is designed to be run inside the AoU Researcher Workbench


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in cromshell cromshell/inputs cromshell/inputs/templates \
           cromshell/job_ids cromshell/progress staging staging/misc; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}
find code/ -name "*.R" | xargs -I {} chmod a+x {}

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh

# Install necessary packages
. code/refs/install_packages.sh python R

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Infer workspace number and save as environment variable
export WN=$( get_workspace_number )

# Download workspace-specific contig lists
gsutil cp -r \
  gs://dfci-g2c-refs/hg38/contig_lists \
  ./


############################################
# Collect initial short variant QC metrics #
############################################

# Write template .json of inputs for chromsharded manager
# TODO: implement this

# Submit, monitor, stage, and cleanup short variant QC metadata task
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/CollectVcfQcMetrics.wdl \
  --input-json-template TBD \
  --contig-variable-overrides TBD \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/ \
  --name CollectShortVariantQcMetrics \
  --contig-list contig_lists/dfci-g2c.v1.contigs.w$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.initial_qc.CollectShortVariantQcMetrics.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1.gatkhc" \
  --outer-gate 30 \
  --max-attempts 2


#################################
# Collect initial SV QC metrics #
#################################

# Write template .json of inputs for chromsharded manager
# TODO: implement this

# Submit, monitor, stage, and cleanup SV QC metadata task
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/CollectVcfQcMetrics.wdl \
  --input-json-template TBD \
  --contig-variable-overrides TBD \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/ \
  --name CollectSvQcMetrics \
  --contig-list contig_lists/dfci-g2c.v1.contigs.w$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.initial_qc.CollectSvQcMetrics.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1.gatksv" \
  --outer-gate 30 \
  --max-attempts 2


##########################################
# Analyze & visualize initial QC metrics #
##########################################

# TODO: implement this


