#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to perform intake QC & batching for DFCI G2C phase 1 cohort

# Note that this code is designed to be run inside the AoU Researcher Workbench


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in data data/tarballs plots plots/raw_intake_qc; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
for suffix in py R; do
  find code/ -name "*.${suffix}" | xargs -I {} chmod a+x {}
done

# Source .bashrc and bash utility functions
. ~/code/refs/dotfiles/aou.rw.bashrc
. code/refs/aou_bash_utils.sh

# Install necessary packages
. code/refs/install_packages.sh R

# Localize intake QC data and references
gsutil -m -u $GPROJECT cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/dfci-g2c.intake_qc.aou.tsv.gz \
  gs://dfci-g2c-inputs/intake_qc/dfci-g2c.intake_qc.non_aou.tsv.gz \
  gs://dfci-g2c-inputs/phenotypes/dfci-g2c.non_aou.phenos.tsv.gz \
  gs://dfci-g2c-refs/hgsv/HGSV.ByrskaBishop.sample_populations.tsv \
  data/
gsutil -m -u $GPROJECT cat \
  $MAIN_WORKSPACE_BUCKET/data/sample_info/other/all_g2c.phase1.samples.list \
  gs://dfci-g2c-inputs/intake_qc/all_g2c.phase1.samples.non_aou.list \
| sort -V > data/all_g2c.phase1.intake_qc.samples.list


########################
# INTAKE QC & BATCHING #
########################

# Add minimal necessary phenotype data to non-AoU and AoU samples for QC + batching
for subset in non_aou aou; do
  code/scripts/join_qc_and_phenotypes.R \
    --qc-tsv data/dfci-g2c.intake_qc.$subset.tsv.gz \
    --phenos-tsv data/dfci-g2c.$subset.phenos.tsv.gz \
    --out-tsv data/dfci-g2c.intake_qc_with_phenos.$subset.tsv
  gzip -f data/dfci-g2c.intake_qc_with_phenos.$subset.tsv
done

# Merge AoU and non-AoU manifests
code/scripts/merge_intake_qc_manifests.R \
  --aou-tsv data/dfci-g2c.intake_qc_with_phenos.aou.tsv.gz \
  --non-aou-tsv data/dfci-g2c.intake_qc_with_phenos.non_aou.tsv.gz \
  --include-samples data/all_g2c.phase1.intake_qc.samples.list \
  --hgsv-pop-assignments data/HGSV.ByrskaBishop.sample_populations.tsv \
  --id-prefix G2C \
  --suffix-length 6 \
  --out-tsv data/dfci-g2c.intake_qc.all.tsv
gzip -f data/dfci-g2c.intake_qc.all.tsv

# Copy merged manifest to main workspace bucket
gsutil -m cp \
  data/dfci-g2c.intake_qc.all.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/

# Generate initial QC plots prior to batching
code/scripts/plot_intake_qc.R \
  --qc-tsv data/dfci-g2c.intake_qc.all.tsv.gz \
  --out-prefix plots/raw_intake_qc/dfci-g2c.phase1.raw_intake_qc

# Compress & copy initial QC plots to results bucket
tar -czvf data/tarballs/dfci-g2c.phase1.raw_intake_qc.plots.tar.gz plots/raw_intake_qc
gsutil -m cp \
  data/tarballs/dfci-g2c.phase1.raw_intake_qc.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/intake_qc/

# Run batching & QC procedure
# TODO: implement this

# Generate QC plots after batching (on full cohort, not per-batch)
# TODO: implement this

# Generate summary plots per batch
# TODO: implement this

