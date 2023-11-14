#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Cleans up a single batch of HMF germline WGS data processed in GCP

# Note: many of the buckets listed below have strict access permissions
# This script will not work unless you have read/write permissions to all buckets

MDAT=$1
BATCH=$2

usage() {
cat << EOF

USAGE: ./cleanup_hmf_batch.sh sample_info.tsv batch_name

       Handles all data organization & cloud storage cleaning for a single
       batch of HMF samples processed in Terra/GCP

       Note that extensive GCP permissions are required for this script to function

EOF
}

if [ $# -lt 2 ]; then
  echo -e "\nERROR: must supply Terra-style sample metadata .tsv and batch name as positional arguments"
  usage
  exit
fi

set -eu -o pipefail

export ARCHIVAL_BUCKET="gs://dfci-g2c-archival"
export STAGING_BUCKET="gs://dfci-g2c-inputs"
export SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )


# Set up files
WRKDIR=`mktemp -d`
cat \
  <( head -n1 $MDAT ) \
  <( fgrep -w $BATCH $MDAT ) \
> $WRKDIR/mdat.tsv
sed '1d' $WRKDIR/mdat.tsv | cut -f1 > $WRKDIR/samples.list
for colname in $( head -n1 $WRKDIR/mdat.tsv | cut -f2- ); do
  cidx=$( head -n1 $WRKDIR/mdat.tsv | sed 's/\t/\n/g' | awk -v cname=$colname '{ if ($1==cname) print NR }' )
  sed '1d' $WRKDIR/mdat.tsv | cut -f$cidx > $WRKDIR/$colname.list
done


######################################################
# 1. Relocate all outputs to G2C inputs storage bucket
######################################################

# hg38-aligned CRAMs
for field in bam_or_cram_file bam_or_cram_index; do
  cat $WRKDIR/$field.list
done | gsutil -m mv -I $ARCHIVAL_BUCKET/crams/hmf/

# GATK-HC gVCFs
for field in old_gvcf old_gvcf_index; do
  cat $WRKDIR/$field.list
done | gsutil -m mv -I $STAGING_BUCKET/hmf/gatk-hc/
for field in reblocked_gvcf reblocked_gvcf_idx; do
  cat $WRKDIR/$field.list
done | gsutil -m mv -I $STAGING_BUCKET/hmf/gatk-hc/reblocked/

# Manta
for field in manta_vcf manta_index; do
  cat $WRKDIR/$field.list
done | gsutil -m mv -I $STAGING_BUCKET/hmf/manta/

# MELT
for field in melt_vcf melt_index; do
  cat $WRKDIR/$field.list
done | gsutil -m mv -I $STAGING_BUCKET/hmf/melt/

# Wham
for field in wham_vcf wham_index; do
  cat $WRKDIR/$field.list
done | gsutil -m mv -I $STAGING_BUCKET/hmf/wham/

# GATK-SV files
cat $WRKDIR/coverage_counts.list \
| gsutil -m mv -I $STAGING_BUCKET/hmf/gatk-sv/coverage/
for field in pesr_disc pesr_disc_index pesr_split pesr_split_index pesr_sd pesr_sd_index; do
  cat $WRKDIR/$field.list
done | gsutil -m mv -I $STAGING_BUCKET/hmf/gatk-sv/pesr/
cat $WRKDIR/sample_metrics_files.list | sed 's/,/\n/g' | tr -d '"[]' \
| gsutil -m mv -I $STAGING_BUCKET/hmf/gatk-sv/metrics/


#####################################
# 2. Check for completion of transfer
#####################################
gsutil -m ls $STAGING_BUCKET/hmf/** > $WRKDIR/all_uris.list
while read sample; do
  $SCRIPT_DIR/../sample_info/check_sample_inputs.sh \
    hmf $sample $WRKDIR/all_uris.list
done < $WRKDIR/samples.list > $WRKDIR/sample_completion.txt
n_incomplete=$( awk '{ if ($2!="") print }' $WRKDIR/sample_completion.txt | wc -l )
if [ $n_incomplete -gt 0 ]; then
  echo -e "\n\nERROR: $n_incomplete samples incomplete. Exiting."
  cat $WRKDIR/sample_completion
  exit 1
fi


##############################################################################
# 3. Clear all temporary processing execution buckets for all samples in batch
##############################################################################
# USA-hosted hg19 CRAMs
for field in hg19_cram_usa hg19_crai_usa; do
  cat $WRKDIR/$field.list
done | gsutil -m rm -I

# Terra-generated files
# Note: focus on large files (not logs etc) to make cleaning faster
cat $WRKDIR/mdat.tsv | sed 's/[\t,]/\n/g' | fgrep "gs://fc-secure-" \
| tr -d '"[]' | cut -d/ -f1-7 | sort -V | uniq | awk '{ print $1"/**" }' \
> $WRKDIR/gsdirs_to_rm.list
gsutil -m ls $( paste -s -d\  $WRKDIR/gsdirs_to_rm.list ) \
| grep '\.gz$\|\.bam$\|\.cram$|\.txt$|\.tsv$|\.bgz$' \
> $WRKDIR/uris_to_rm.list || true
cat $WRKDIR/uris_to_rm.list | gsutil -m rm -I


############
# 4. Local cleanup
############
rm -rf $WRKDIR

