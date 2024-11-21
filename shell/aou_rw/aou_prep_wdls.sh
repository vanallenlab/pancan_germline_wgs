#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Code to copy all necessary WDLs and other Cromwell-related files (e.g., 
# input .json templates) to AoU RW bucket (as RW has no internet access)

# See aou_prep_libs.sh for libraries & tools to be executed on RW terminal directly

# Note that this code is designed to be run *locally* (not on RW)

# Set up local working directory
WRKDIR=`mktemp -d`
cd $WRKDIR

# Clone G2C repo & checkout branch of interest
export g2c_branch=gatksv
git clone git@github.com:vanallenlab/pancan_germline_wgs.git --branch=$g2c_branch

# Clone GATK-SV repo & checkout release tag of interest
export gatksv_tag=v1.0
git clone git@github.com:broadinstitute/gatk-sv.git --branch=$gatksv_tag

# Clone GATK-HC workflows repo & checkout release tag of interest
export gatkhc_tag=2.3.1
git clone git@github.com:gatk-workflows/gatk4-germline-snps-indels.git --branch=$gatkhc_tag

# Make & populate directory of all WDLs and other reference files
for dir in wdl wdl/pancan_germline_wgs wdl/gatk-sv wdl/gatk-hc; do
  if [ -e $dir ]; then rm -rf $dir; fi
  mkdir $dir
done
cp pancan_germline_wgs/wdl/*.wdl $WRKDIR/wdl/pancan_germline_wgs/
cp gatk-sv/wdl/*.wdl $WRKDIR/wdl/gatk-sv/
cp gatk4-germline-snps-indels/*.wdl $WRKDIR/wdl/gatk-hc/

# Copy WDLs to AoU RW bucket
# Note: must use AoU Google credentials
export rw_bucket=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45
# gcloud auth login
gsutil -m cp -r \
  wdl \
  pancan_germline_wgs/refs \
  $rw_bucket/code/

# Clean up
cd ~
rm -rf $WRKDIR
