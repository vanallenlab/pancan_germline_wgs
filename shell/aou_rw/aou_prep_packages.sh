#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Code to copy all necessary analysis scripts and software packages
# to AoU RW bucket (as RW has no internet access)

# Note that this code is designed to be run *locally* (not on RW)

# Set up local working directory
WRKDIR=`mktemp -d`
cd $WRKDIR

# Clone G2C repo & checkout branch of interest
export g2c_branch=aou_processing
git clone git@github.com:vanallenlab/pancan_germline_wgs.git && \
cd pancan_germline_wgs && \
git checkout $g2c_branch && \
git pull && \
cd ../

# Clone RLCtools repo & checkout branch of interest
export rlctools_branch=main
git clone git@github.com:RCollins13/RLCtools.git && \
cd RLCtools && \
git checkout $rlctools_branch && \
git pull && \
cd ../

# Copy code to AoU RW bucket
# Note: must use AoU Google credentials
export rw_bucket=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45
gcloud auth login
gsutil -m cp -r \
  pancan_germline_wgs/src/G2C_*.tar.gz \
  RLCtools/RLCtools_*.tar.gz \
  $rw_bucket/code/src/

# Clean up
cd ~
rm -rf $WRKDIR
