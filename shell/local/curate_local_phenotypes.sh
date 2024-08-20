#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to curate selected sample phenotypes for DFCI G2C phase 1 cohort

# Note that this code is designed to be run on RLC's local machine
# and references local paths that will not generalize to other systems

# This code is intended for record-keeping/reproducibility purposes
# rather than reapplication by other users


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export WRKDIR=`mktemp -d`
export BASEDIR="/Users/ryan/Desktop/Collins/VanAllen/pancancer_wgs"
export CODEDIR="$BASEDIR/pancan_germline_wgs"


###########
# NCI GDC #
###########

# Apollo
export APOLLODIR="$BASEDIR/data_and_cohorts/apollo/clinical.cases_selection.2024-03-26"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $APOLLODIR/clinical.tsv \
  --exposure-tsv $APOLLODIR/exposure.tsv \
  --biospecimen-tsv $APOLLODIR/../biospecimen.cases_selection.2024-03-26/sample.tsv \
  --cohort apollo \
  --out-tsv $WRKDIR/apollo.phenos.tsv

# CPTAC
export CPTACDIR="$BASEDIR/data_and_cohorts/cptac/clinical.cases_selection.2024-03-27"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $CPTACDIR/clinical.tsv \
  --exposure-tsv $CPTACDIR/exposure.tsv \
  --biospecimen-tsv $CPTACDIR/../biospecimen.project-cptac-3.2024-08-20/sample.tsv \
  --cohort cptac \
  --out-tsv $WRKDIR/cptac.phenos.tsv

# Eagle
export EAGLEDIR="$BASEDIR/data_and_cohorts/eagle/clinical.cases_selection.2024-03-29"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $EAGLEDIR/clinical.tsv \
  --exposure-tsv $EAGLEDIR/exposure.tsv \
  --biospecimen-tsv $EAGLEDIR/../biospecimen.cases_selection.2024-03-29/sample.tsv \
  --cohort eagle \
  --out-tsv $WRKDIR/eagle.phenos.tsv

# HCMI
export HCMIDIR="$BASEDIR/data_and_cohorts/hcmi/clinical.cases_selection.2024-03-26"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $HCMIDIR/clinical.tsv \
  --exposure-tsv $HCMIDIR/exposure.tsv \
  --biospecimen-tsv $HCMIDIR/../biospecimen.cases_selection.2024-03-26/sample.tsv \
  --cohort hcmi \
  --out-tsv $WRKDIR/hcmi.phenos.tsv

# WCDT
export WCDTDIR="$BASEDIR/data_and_cohorts/wcdt/clinical.cases_selection.2024-03-26"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $WCDTDIR/clinical.tsv \
  --exposure-tsv $WCDTDIR/exposure.tsv \
  --biospecimen-tsv $WCDTDIR/..//biospecimen.cases_selection.2024-03-26/sample.tsv \
  --cohort wcdt \
  --out-tsv $WRKDIR/wcdt.phenos.tsv


#########
# dbGaP #
#########

# BioMe
# TODO: add this

# CEPH
# TODO: add this

# GTEx
# TODO: add this

# MESA
# TODO: add this

# LCINS
# TODO: add this


#################
# OTHER|BESPOKE #
#################

# HGSVC
# TODO: add this

# HMF
# TODO: add this

# ICGC
# TODO: add this

# Proactive
# TODO: add this

# UFC
# TODO: add this


###################
# STAGE & CLEANUP #
###################

# Pool all phenotype data
head -n1 $( find $WRKDIR/ -name "*.phenos.tsv" | head -n1 ) \
> $WRKDIR/dfci-g2c.non_aou.phenos.tsv
find $WRKDIR/ -name "*.phenos.tsv" \
| xargs -I {} sed '1d' {} \
| sort -Vk2,2 -k1,1V \
>> $WRKDIR/dfci-g2c.non_aou.phenos.tsv

# Compress all phenotype data
find $WRKDIR/ -name "*.phenos.tsv" | xargs -I {} gzip -f {}

# Copy all per-cohort phenotype files to main (non-AoU) project bucket for storage
gsutil -m cp \
  $WRKDIR/*.phenos.tsv.gz \
  gs://dfci-g2c-inputs/phenotypes/

# Clean up working directory
rm -rf $WRKDIR

