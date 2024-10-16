#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Check Google cloud storage bucket for availability of all data for a single sample ID

COHORT=$1
SAMPLE=$2
BUCKET=$3  # Optional base bucket to query
GZ_LIST=$4 # Optional pre-computed list of GCP URIs to query

usage() {
cat << EOF

USAGE: ./check_sample_inputs.sh cohort sample_id [bucket] [uri_list]

       Checks google cloud storage bucket for complete input data for a single sample ID

       Can optionally provide:
       - Custom bucket URI
       - Pre-computed list of GCP URIs, which will be referenced directly 
         instead of querying bucket on-the-fly

EOF
}

if [ $# -lt 2 ]; then
  echo -e "\nERROR: must supply cohort and sample name as positional arguments"
  usage
  exit
fi

if [ -z $BUCKET ]; then
  BUCKET="gs://dfci-g2c-inputs"
fi

if [ -z $TMPDIR ]; then
  TMPDIR=`mktemp -d`
  RM_TMPDIR="true"
else
  RM_TMPDIR="false"
fi

# Build list of files that are expected to exist for each sample
export TAB=$'\t'
cat << EOF > $TMPDIR/$COHORT.$SAMPLE.paths
gatksv_coverage$TAB$BUCKET/$COHORT/gatk-sv/coverage/$SAMPLE.counts.tsv.gz
gatksv_pe$TAB$BUCKET/$COHORT/gatk-sv/pesr/$SAMPLE.pe.txt.gz
gatksv_sd$TAB$BUCKET/$COHORT/gatk-sv/pesr/$SAMPLE.sd.txt.gz
gatksv_sr$TAB$BUCKET/$COHORT/gatk-sv/pesr/$SAMPLE.sr.txt.gz
gatkhc_reblocked_gvcf$TAB$BUCKET/$COHORT/gatk-hc/reblocked/$SAMPLE.reblocked.g.vcf.gz
manta$TAB$BUCKET/$COHORT/manta/$SAMPLE.manta.vcf.gz
melt$TAB$BUCKET/$COHORT/melt/$SAMPLE.melt.vcf.gz
wham$TAB$BUCKET/$COHORT/wham/$SAMPLE.wham.vcf.gz
demographics$TAB$BUCKET/$COHORT/demographics/$SAMPLE.txt
charr$TAB$BUCKET/$COHORT/charr/$SAMPLE.txt
EOF

# Check which files are present in gs://
if [ -z $GZ_LIST ]; then
  gsutil -m ls $( awk -v ORS=" " '{ print $2 }' $TMPDIR/$COHORT.$SAMPLE.paths ) \
  1> $TMPDIR/$COHORT.$SAMPLE.objects_found 2> /dev/null
else
  cp $GZ_LIST $TMPDIR/$COHORT.$SAMPLE.objects_found
fi

# Build report of missing files
missing=""
while read dname path; do
  if [ $( awk -v path="$path" '{ if ($1==path) print }' \
          $TMPDIR/$COHORT.$SAMPLE.objects_found | wc -l ) -eq 0 ]; then
    missing="${missing},$dname"
  fi
done < $TMPDIR/$COHORT.$SAMPLE.paths
missing=$( echo $missing | sed 's/,/\n/g' | sed '/^$/d' | sort -V | paste -s -d, )
echo -e "$SAMPLE\t$missing"

# Clean up
rm $TMPDIR/$COHORT.$SAMPLE.paths $TMPDIR/$COHORT.$SAMPLE.objects_found
if [ $RM_TMPDIR == "true" ]; then rm -rf $TMPDIR; fi
  