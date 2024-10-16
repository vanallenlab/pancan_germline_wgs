#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins, Noah Fields, and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This WDL is the refactored successor of "Consolidate_Cohort_Data.wdl", originally by Noah Fields
# It collects intake sample QC information into a single table for a single cohort

COHORT=$1
PLOIDY=$2
BUCKET=$3  # Optional base bucket to query

usage() {
cat << EOF

USAGE: ./consolidate_cohort_intake_qc.sh cohort ploidy.tsv [bucket]

       Collects intake sample QC information into a single table for a single cohort

       Can optionally provide custom bucket URI as third positional argument

EOF
}

# Parse positional arguments
if [ $# -lt 2 ]; then
  echo -e "\nERROR: must supply cohort name and ploidy.tsv as positional arguments"
  usage
  exit
fi

if [ -z $BUCKET ]; then
  BUCKET="gs://dfci-g2c-inputs"
fi

# Ensure $TMPDIR is set (files will be staged here)
if [ -z $TMPDIR ]; then
  TMPDIR=`mktemp -d`
  RM_TMPDIR="true"
else
  RM_TMPDIR="false"
fi

# Concatenate demographic data
mkdir $TMPDIR/${COHORT}_demo
gsutil -m cp ${BUCKET}/${COHORT}/demographics/*.txt $TMPDIR/${COHORT}_demo/
head -n1 $( find $TMPDIR/${COHORT}_demo/ -name "*.txt" | head -n1 ) \
> $TMPDIR/${COHORT}_demo/demo.header
cat $TMPDIR/${COHORT}_demo/*.txt | grep -v '^Sample' | sed '/^$/d' \
| sort -Vk1,1 | uniq | cat $TMPDIR/${COHORT}_demo/demo.header - | gzip -c \
> $TMPDIR/${COHORT}_demo/$COHORT.demo.tsv.gz

# Concatenate Charr data
mkdir $TMPDIR/${COHORT}_charr
gsutil -m cp ${BUCKET}/${COHORT}/charr/*.txt $TMPDIR/${COHORT}_charr/
head -n1 $( find $TMPDIR/${COHORT}_charr/ -name "*.txt" | head -n1 ) \
> $TMPDIR/${COHORT}_charr/charr.header
cat $TMPDIR/${COHORT}_charr/*.txt | grep -v '^#SAMPLE' | sed '/^$/d' \
| sort -Vk1,1 | uniq | cat $TMPDIR/${COHORT}_charr/charr.header - | gzip -c \
> $TMPDIR/${COHORT}_charr/$COHORT.charr.tsv.gz

# Concatenate GATK-SV coverage metrics
mkdir $TMPDIR/${COHORT}_cov
gsutil -m cp \
  ${BUCKET}/${COHORT}/gatk-sv/metrics/*.raw-counts.tsv \
  $TMPDIR/${COHORT}_cov/
find $TMPDIR/${COHORT}_cov/ -name "*.tsv" > $TMPDIR/${COHORT}_cov/cov.files.list
paste \
  <( cat $TMPDIR/${COHORT}_cov/cov.files.list | xargs -I {} basename {} | sed 's/.raw-counts.tsv//g' ) \
  <( cat $TMPDIR/${COHORT}_cov/cov.files.list | xargs -I {} cat {} | grep "^rd_q50_" | cut -f2 ) \
  <( cat $TMPDIR/${COHORT}_cov/cov.files.list | xargs -I {} cat {} | grep "^rd_mean_" | cut -f2 ) \
| sort -Vk1,1 | cat <( echo -e "Sample\trd_median\trd_mean" ) - | gzip -c \
> $TMPDIR/${COHORT}_cov/$COHORT.cov.tsv.gz

# Copy ploidy data
gsutil -m cp $PLOIDY $TMPDIR/$COHORT.ploidy.tsv.gz

# Merge all data, keeping only the strict intersection of sample IDs present in
# Charr, demographics, and ploidy. It will tolerate missingness in GATK-SV coverage.
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
$SCRIPT_DIR/merge_cohort_intake_qc_components.R \
  --charr $TMPDIR/${COHORT}_charr/$COHORT.charr.tsv.gz \
  --demo $TMPDIR/${COHORT}_demo/$COHORT.demo.tsv.gz \
  --cov $TMPDIR/${COHORT}_cov/$COHORT.cov.tsv.gz \
  --ploidy $TMPDIR/$COHORT.ploidy.tsv.gz \
  --outfile $TMPDIR/$COHORT.intake_qc.tsv
gzip -f $TMPDIR/$COHORT.intake_qc.tsv

# Copy merged intake QC metrics to $BUCKET
gsutil -m cp \
  $TMPDIR/$COHORT.intake_qc.tsv.gz \
  $BUCKET/intake_qc/cohort_metrics/

# Clean up
rm -rf \
  $TMPDIR/${COHORT}_charr \
  $TMPDIR/${COHORT}_demo \
  $TMPDIR/${COHORT}_cov \
  $TMPDIR/$COHORT.intake_qc.tsv.gz
if [ $RM_TMPDIR == "true" ]; then rm -rf $TMPDIR; fi
