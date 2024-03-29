#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Check Google cloud storage bucket for availability of all data for all cohorts

set -eu -o pipefail

BUCKET="gs://dfci-g2c-inputs"
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
WRKDIR=`mktemp -d`

usage() {
cat << EOF

USAGE: ./check_all_sample_inputs.sh [-h|--help] [COHORTS]

       Checks google cloud storage bucket for complete input data for all DFCI-G2C cohorts

       Options:
       COHORTS: comma-delimited list of cohorts to check [default: check all cohorts]

EOF
}

if [ $# -lt 1 ]; then
  gsutil ls $BUCKET/sample-lists/*.samples.list \
  | xargs -I {} basename {} \
  | sed 's/\.samples\.list//g' \
  > $WRKDIR/cohorts.list
elif [ $1 == "--help" ] || [ $1 == "-h" ]; then
  usage
  exit 0
else
  echo $1 | sed 's/,/\n/g' | sort -V | uniq > $WRKDIR/cohorts.list
fi

k=0
while read COHORT; do
  nsamp=$( gsutil cat $BUCKET/sample-lists/$COHORT.samples.list | wc -l )
  echo -e "Checking status of $nsamp samples for $COHORT"

  # For efficiency, first compile list of all URIs suffixed with *.gz
  gsutil -m ls $BUCKET/$COHORT/**.gz \
  1> $WRKDIR/$COHORT.objects_found 2> /dev/null

  echo -e "#sample\tmissing_inputs" > $WRKDIR/$COHORT.sample_status.tsv
  while read SAMPLE; do
    k=$((k+1))
    echo -e "Sample $k: $SAMPLE"
    $SCRIPT_DIR/check_sample_inputs.sh \
      $COHORT $SAMPLE $WRKDIR/$COHORT.objects_found \
    >> $WRKDIR/$COHORT.sample_status.tsv
  done < <( gsutil cat $BUCKET/sample-lists/$COHORT.samples.list )
done < <( sed '/^$/d' $WRKDIR/cohorts.list )

# Localize status to bucket and remove local copy
gsutil cp $WRKDIR/*.tsv $BUCKET/sample-status/
rm -rf $WRKDIR
