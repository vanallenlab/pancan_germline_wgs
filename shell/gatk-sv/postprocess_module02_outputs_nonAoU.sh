#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Local helper code to process GATK-SV module 02 (ploidy/WGD) for all non-AoU samples

# Note that this code is designed to be run locally on RLC's computer 
# and references paths that won't generalize to other systems


# Set paths
export WRKDIR=/Users/ryan/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/ploidy_wgd
export CODEDIR=/Users/ryan/Desktop/Collins/VanAllen/pancancer_wgs/pancan_germline_wgs/
export TMPDIR=`mktemp -d`
cd $WRKDIR


# First, get indexes for relevant columns from manifest
ploidy_plots_idx=$( head -n1 $WRKDIR/dfci-g2c-ploidy-estimation.terra_manifest.tsv \
                    | sed 's/\t/\n/g' \
                    | awk -v query="ploidy_plots" '{ if ($1==query) print NR }' )
wgd_plot_idx=$( head -n1 $WRKDIR/dfci-g2c-ploidy-estimation.terra_manifest.tsv \
                | sed 's/\t/\n/g' \
                | awk -v query="WGD_dist" '{ if ($1==query) print NR }' )
qc_table_idx=$( head -n1 $WRKDIR/dfci-g2c-ploidy-estimation.terra_manifest.tsv \
                | sed 's/\t/\n/g' \
                | awk -v query="qc_table" '{ if ($1==query) print NR }' )
wgd_tsv_idx=$( head -n1 $WRKDIR/dfci-g2c-ploidy-estimation.terra_manifest.tsv \
               | sed 's/\t/\n/g' \
               | awk -v query="WGD_scores" '{ if ($1==query) print NR }' )
med_cov_idx=$( head -n1 $WRKDIR/dfci-g2c-ploidy-estimation.terra_manifest.tsv \
               | sed 's/\t/\n/g' \
               | awk -v query="bincov_median" '{ if ($1==query) print NR }' )

# Process each batch in serial
while read bid; do

  # Make directory
  if ! [ -e $WRKDIR/$bid ]; then
    mkdir $WRKDIR/$bid
  fi

  # For convenience, write list of samples/cohorts from batch to batch directory
  awk -v bid=$bid '{ if ($1==bid) print $2 }' \
    $WRKDIR/dfci-g2c-ploidy-estimation.batch_membership.tsv \
  | sed 's/_/\t/' | sort -Vk1,1 -k2,2V \
  | awk -v OFS="\t" '{ print $2, $1 }' \
  | cat <( echo -e "#sample\tcohort" ) - \
  > $WRKDIR/$bid/$bid.sample_info.tsv

  # Download median coverage, WGD plot, and WGD .tsv for permanent local storage
  awk -v bid=$bid -v idx1=$wgd_plot_idx -v idx2=$wgd_tsv_idx -v idx3=$med_cov_idx \
    -v OFS="\n" -v FS="\t" '{ if ($1==bid) print $idx1, $idx2, $idx3 }' \
    $WRKDIR/dfci-g2c-ploidy-estimation.terra_manifest.tsv \
  | gsutil -m cp -I $WRKDIR/$bid/

  # Download & clean QC table
  awk -v bid=$bid -v idx1=$qc_table_idx -v OFS="\n" -v FS="\t" \
    '{ if ($1==bid) print $idx1 }' \
    $WRKDIR/dfci-g2c-ploidy-estimation.terra_manifest.tsv \
  | gsutil -m cp -I $TMPDIR/
  $CODEDIR/scripts/sample_info/clean_module02_qc_table.R \
    $TMPDIR/$bid.evidence_qc_table.tsv
  mv $TMPDIR/$bid.evidence_qc_table.tsv $WRKDIR/$bid/

  # Download & process ploidy plots tarball (don't need to save everything)
  pp_uri=$( awk -v bid=$bid -v idx1=$ploidy_plots_idx -v OFS="\n" -v FS="\t" \
              '{ if ($1==bid) print $idx1 }' \
              $WRKDIR/dfci-g2c-ploidy-estimation.terra_manifest.tsv )
  gsutil -m cp $pp_uri $TMPDIR/
  tar -xzvf $TMPDIR/$( basename $pp_uri ) --directory $TMPDIR/
  find $TMPDIR/ploidy_est/*png | xargs -I {} mv {} $WRKDIR/$bid/
  mv $TMPDIR/ploidy_est/sample_sex_assignments.txt.gz $WRKDIR/$bid/
  rm -rf $TMPDIR/ploidy_est $TMPDIR/$( basename $pp_uri )

done < <( sed '1d' $WRKDIR/dfci-g2c-ploidy-estimation.terra_manifest.tsv | cut -f1 )

# Once all batches have been processed, combine the QC tables across all batches
head -n1 $( find $WRKDIR -name "*.evidence_qc_table.tsv" | head -n1 ) \
| sed 's/^#ID/#cohort\tsample/' \
> $WRKDIR/dfci-g2c-ploidy-estimation.merged_qc_table.tsv
find $WRKDIR -name "*.evidence_qc_table.tsv" \
| xargs -I {} cat {} | grep -ve '^#' | sed 's/_/\t/' | sort -Vk1,1 -k2,2V \
>> $WRKDIR/dfci-g2c-ploidy-estimation.merged_qc_table.tsv
gzip -f $WRKDIR/dfci-g2c-ploidy-estimation.merged_qc_table.tsv

# Copy QC table to G2C staging bucket
gsutil -m cp \
  $WRKDIR/dfci-g2c-ploidy-estimation.merged_qc_table.tsv.gz \
  gs://dfci-g2c-inputs/intake_qc/


# Cleanup
rm -rf $TMPDIR
