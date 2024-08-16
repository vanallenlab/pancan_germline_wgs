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
export MAIN_AOU_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45
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
  --cohort apollo \
  --out-tsv $WRKDIR/apollo.phenos.tsv

# CPTAC
export CPTACDIR="$BASEDIR/data_and_cohorts/cptac/clinical.cases_selection.2024-03-27"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $CPTACDIR/clinical.tsv \
  --exposure-tsv $CPTACDIR/exposure.tsv \
  --cohort cptac \
  --out-tsv $WRKDIR/cptac.phenos.tsv

# Eagle
export EAGLEDIR="$BASEDIR/data_and_cohorts/eagle/clinical.cases_selection.2024-03-29"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $EAGLEDIR/clinical.tsv \
  --exposure-tsv $EAGLEDIR/exposure.tsv \
  --cohort eagle \
  --out-tsv $WRKDIR/eagle.phenos.tsv

# HCMI
export HCMIDIR="$BASEDIR/data_and_cohorts/hcmi/clinical.cases_selection.2024-03-26"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $HCMIDIR/clinical.tsv \
  --exposure-tsv $HCMIDIR/exposure.tsv \
  --cohort hcmi \
  --out-tsv $WRKDIR/hcmi.phenos.tsv

# WCDT
export WCDTDIR="$BASEDIR/data_and_cohorts/wcdt/clinical.cases_selection.2024-03-26"
$CODEDIR/scripts/sample_info/phenotypes/curate_nci_gdc_phenotypes.R \
  --clinical-tsv $WCDTDIR/clinical.tsv \
  --exposure-tsv $WCDTDIR/exposure.tsv \
  --cohort wcdt \
  --out-tsv $WRKDIR/wcdt.phenos.tsv

