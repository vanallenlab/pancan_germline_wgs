#!/usr/bin/env bash

# The Genomic Architecture of Human Cancers
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Check Google cloud storage bucket for availability of all data for a single sample ID

COHORT=$1
SAMPLE=$2

usage() {
cat << EOF

USAGE: ./check_sample_inputs.sh cohort sample_id

       Checks google cloud storage bucket for complete input data for a single sample ID

EOF
}

if [ $# -ne 2 ]; then
  echo -e "\nERROR: must supply cohort and sample name as positional arguments"
  usage
  exit
fi

BUCKET="gs://dfci-g2c-inputs"

# Build list of files that are expected to exist for each sample
export TAB=$'\t'
cat << EOF > $TMPDIR/$COHORT.$SAMPLE.paths
gatksv_coverage$TAB$BUCKET/$COHORT/gatk-sv/coverage/$SAMPLE.counts.tsv.gz
gatksv_pe$TAB$BUCKET/$COHORT/gatk-sv/pesr/$SAMPLE.pe.txt.gz
gatksv_sd$TAB$BUCKET/$COHORT/gatk-sv/pesr/$SAMPLE.sd.txt.gz
gatksv_sr$TAB$BUCKET/$COHORT/gatk-sv/pesr/$SAMPLE.sr.txt.gz
gatkhc_gvcf$TAB$BUCKET/$COHORT/gatk-hc/$SAMPLE.g.vcf.gz
manta$TAB$BUCKET/$COHORT/manta/$SAMPLE.manta.vcf.gz
melt$TAB$BUCKET/$COHORT/melt/$SAMPLE.melt.vcf.gz
wham$TAB$BUCKET/$COHORT/wham/$SAMPLE.wham.vcf.gz
EOF

# Check which files are present in gs://
gsutil -m ls $( awk -v ORS=" " '{ print $2 }' $TMPDIR/$COHORT.$SAMPLE.paths ) \
1> $TMPDIR/$COHORT.$SAMPLE.objects_found 2> /dev/null

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