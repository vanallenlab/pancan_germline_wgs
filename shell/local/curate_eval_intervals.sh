#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to curate QC evaluation intervals from Genome in a Bottle

# Note that this code is designed to be run on RLC's local machine
# and references local paths that will not generalize to other systems

# This code is intended for record-keeping/reproducibility purposes
# rather than reapplication by other users


# Set up working directory
export WRKDIR=`mktemp -d`

# Download Broad hg38 calling intervals and convert to BED
gsutil -m cat \
  gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list \
| grep -ve "^@" \
| cut -f1-3 \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bgzip -c \
> $WRKDIR/broad.hg38.gatkhc_callable.bed.gz

# Invert calling intervals to produce non-callable hg38 intervals
gsutil -m cat gs://dfci-g2c-refs/hg38/hg38.genome \
| awk -v OFS="\t" '{ print $1, "0", $2 }' \
| bedtools subtract -a - -b $WRKDIR/broad.hg38.gatkhc_callable.bed.gz \
| bgzip -c \
> $WRKDIR/broad.hg38.gatkhc_not_callable.bed.gz

# Download GIAB hg38 difficult intervals and keep the intersection vs. GATK-HC callable regions
wget -O - \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/Union/GRCh38_alldifficultregions.bed.gz \
| gunzip -c \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bedtools subtract -a - -b $WRKDIR/broad.hg38.gatkhc_not_callable.bed.gz \
| bgzip -c \
> $WRKDIR/giab.hg38.broad_callable.hard.bed.gz
tabix -p bed -f $WRKDIR/giab.hg38.broad_callable.hard.bed.gz

# Download GIAB hg38 easy intervals and keep the intersection vs. GATK-HC callable regions
wget -O - \
  https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/genome-stratifications/v3.5/GRCh38@all/Union/GRCh38_notinalldifficultregions.bed.gz \
| gunzip -c \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bedtools subtract -a - -b $WRKDIR/broad.hg38.gatkhc_not_callable.bed.gz \
| bgzip -c \
> $WRKDIR/giab.hg38.broad_callable.easy.bed.gz
tabix -p bed -f $WRKDIR/giab.hg38.broad_callable.easy.bed.gz

# Copy main interval files to GS bucket
gsutil -m cp \
  $WRKDIR/giab.hg38.broad_callable.*.bed.gz* \
  gs://dfci-g2c-refs/giab/

# Split easy & hard eval intervals by chromosome
for k in $( seq 1 22 ) X Y; do
  contig=chr$k
  for suf in easy hard; do
    tabix $WRKDIR/giab.hg38.broad_callable.$suf.bed.gz $contig \
    | bgzip -c \
    > $WRKDIR/giab.hg38.broad_callable.$suf.$contig.bed.gz
    tabix -p bed -f $WRKDIR/giab.hg38.broad_callable.$suf.$contig.bed.gz
  done
  gsutil -m cp \
    $WRKDIR/giab.hg38.broad_callable.*.$contig.bed.gz \
    gs://dfci-g2c-refs/giab/$contig/
done

# Clean up
rm -rf $WRKDIR

