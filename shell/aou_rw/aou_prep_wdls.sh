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
export g2c_branch=posthoc_qc
git clone git@github.com:vanallenlab/pancan_germline_wgs.git --branch=$g2c_branch

# Clone GATK-SV repo & checkout release tag of interest
export gatksv_tag=v1.0.1
git clone git@github.com:broadinstitute/gatk-sv.git --branch=$gatksv_tag

# Clone GATK-HC workflow repos & checkout release tag of interest
export gatkhc_tag=2.3.1
git clone git@github.com:gatk-workflows/gatk4-germline-snps-indels.git --branch=$gatkhc_tag
git clone git@github.com:gatk-workflows/utility-wdls.git
git clone git@github.com:broadinstitute/warp.git

# Make & populate directory of all WDLs and other reference files
for dir in wdl wdl/pancan_germline_wgs wdl/gatk-sv wdl/gatk-hc legacy_mingq_wdl; do
  if [ -e $dir ]; then rm -rf $dir; fi
  mkdir $dir
done
cp pancan_germline_wgs/wdl/*.wdl $WRKDIR/wdl/pancan_germline_wgs/
cp -r pancan_germline_wgs/wdl/vcf-qc $WRKDIR/wdl/pancan_germline_wgs/
cp gatk-sv/wdl/*.wdl $WRKDIR/wdl/gatk-sv/
cp gatk4-germline-snps-indels/*.wdl $WRKDIR/wdl/gatk-hc/
cp utility-wdls/*.wdl $WRKDIR/wdl/gatk-hc/
cp warp/tasks/broad/JointGenotypingTasks.wdl $WRKDIR/wdl/gatk-hc/

# Add legacy version of GATK-SV WDLs from most recent branch with working copy of minGQ
cd gatk-sv && \
git checkout origin/xz_fixes_3_rlc_mod && \
cd - > /dev/null && \
cp gatk-sv/wdl/*.wdl $WRKDIR/legacy_mingq_wdl/

# Override any GATK WDLs with their corresponding custom G2C copies
# This is rarely necessary but was deemed the easiest solution for handling 
# edge cases (like for SV module 06) or situations where we intentionally 
# deviated from GATK default procedures (like for outlier sample definition 
# in SV module 08)
cp pancan_germline_wgs/wdl/gatk-sv/* $WRKDIR/wdl/gatk-sv/
mv $WRKDIR/wdl/gatk-sv/Module07MinGQTasks.wdl $WRKDIR/legacy_mingq_wdl/
cp pancan_germline_wgs/wdl/gatk-hc/* $WRKDIR/wdl/gatk-hc/

# Copy WDLs to AoU RW bucket
# Note: must use AoU Google credentials
export rw_bucket=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45
# gcloud auth login
gsutil -m cp -r \
  wdl \
  pancan_germline_wgs/refs \
  $rw_bucket/code/
gsutil -m cp -r \
  legacy_mingq_wdl \
  $rw_bucket/misc/

# Clean up
cd - >/dev/null
rm -rf $WRKDIR
