#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Check Google cloud storage bucket for availability of all data for all cohorts


usage() {
cat << EOF

USAGE: ./check_all_sample_inputs.sh [COHORT]

       Checks google cloud storage bucket for complete input data for all DFCI-G2C cohorts

       Options:
       COHORT: unless supplied, will check all cohorts

EOF
}

if [ $# -lt 1 ]; then
  cat $TMPDIR/listOfLists.txt | xargs -I {} basename {} | sed 's/\.samples\.list//g' > $TMPDIR/cohorts.list
else
  echo $1 > $TMPDIR/cohorts.list
fi

BUCKET="gs://dfci-g2c-inputs"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

gsutil ls $BUCKET/sample-lists/*.samples.list > $TMPDIR/listOfLists.txt

mkdir $TMPDIR/sample-status
k=0
while read COHORT; do
  nsamp=$( gsutil cat $BUCKET/sample-lists/$COHORT.samples.list | wc -l )
  echo -e "Checking status of $nsamp samples for $COHORT"
  echo -e "#sample\tmissing_inputs" > $TMPDIR/sample-status/$COHORT.sample_status.tsv
  while read SAMPLE; do
    k=$((k+1))
    echo -e "Sample $k: $SAMPLE"
    $SCRIPT_DIR/check_sample_inputs.sh $COHORT $SAMPLE \
    >> $TMPDIR/sample-status/$COHORT.sample_status.tsv
  done < <( gsutil cat $BUCKET/sample-lists/$COHORT.samples.list )
done < $TMPDIR/cohorts.list

# Localize status to bucket and remove local copy
gsutil cp $TMPDIR/sample-status/*.tsv $BUCKET/sample-status/
rm -rf $TMPDIR/sample-status $TMPDIR/cohorts.list