#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Load-balance sets of chromosomes for optimal parallelization across workspaces

# Build local execution directory
WRKDIR=`mktemp -d`
cd $WRKDIR

# TODO: eventually implement this for all five workspaces before scaling up

# For now, we're writing chr21, chr22, and chrY to the "w2" list
# This is only for development purposes of GATK-HC joint genotyping
cat << EOF > dfci-g2c.v1.contigs.w2.list
chr21
chr22
chrY
EOF

# Add three more small contigs to w3 for joint genotyping resource optimization
cat << EOF > dfci-g2c.v1.contigs.w3.list
chr18
chr19
chr20
EOF

# Copy all contig lists to central directory
gsutil -m cp \
  dfci-g2c.v1.contigs.w*.list \
  gs://dfci-g2c-refs/hg38/contig_lists/

# Switch back to execution directory and delete working directory
# (This is only used if/when this script is executed from the command line)
cd -
rm -rf $WRKDIR
