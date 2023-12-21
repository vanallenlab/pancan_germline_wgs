#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-HC on all AoU samples selected for inclusion in G2C pilot

# Note that this code is designed to be run inside the AoU Researcher Workbench

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"

# Prep working directory structure
for dir in data data/cram_paths; do
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

