#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-SV cohort mode pipeline on G2C phase 1

# Note that this code is designed to be run inside the AoU Researcher Workbench
# See gatksv_bash_utils.sh for custom function definitions used below


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in cromshell cromshell/inputs cromshell/job_ids cromshell/progress; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}
find code/ -name "*.R" | xargs -I {} chmod a+x {}

# Install necessary packages
. code/refs/install_packages.sh python

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh
. code/refs/gatksv_bash_utils.sh

# Infer workspace number and save as environment variable
export WN=$( get_workspace_number )

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Create dependencies .zip for all GATK-SV module submissions
cd code/wdl/gatk-sv && \
zip gatksv.dependencies.zip *.wdl && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Copy sample & batch information
gsutil -m cp -r \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/batch_info \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/intake_qc/dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
  ~/


##################
# 03 | TrainGCNV #
##################

module_submission_routine_all_batches 03


############################
# 04 | GatherBatchEvidence #
############################

module_submission_routine_all_batches 04


#####################
# 05 | ClusterBatch #
#####################

module_submission_routine_all_batches 05


##################################
# 05B | ExcludeClusteredOutliers #
##################################

# Note: this is not a canonical GATK-SV module and was instituted specifically for G2C

# First, we need to define a list of outlier samples to be excluded
# This needs to only be run once from one of the AoU workspaces

# Install R packages
. code/refs/install_packages.sh R

# Collect SV count data for each sample per algorithm
if ! [ -e data/gatksv_05B_outliers ]; then
  mkdir gatksv_05B_outliers
fi
for alg in depth manta melt wham; do
  gsutil -m cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/05/*/PlotSVCountsPerSample/CountSVsPerSamplePerType/**.cluster_batch.$alg.svcounts.txt \
  | sort -Vrk1,1 -k2,2V | uniq \
  > gatksv_05B_outliers/dfci-g2c.05_sv_counts.$alg.tsv
done

# Define outlier samples cohort-wide as: 
# * Q3 + 6 IQR for manta, wham, or melt
# * Q3 + 10 IQR for depth CNVs
pop_idx=$( zcat dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
           | sed -n '1p' | sed 's/\t/\n/g' \
           | awk '{ if ($1=="intake_qc_pop") print NR }' )
zcat dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz | sed '1d' \
| awk -v idx=$pop_idx -v FS="\t" -v OFS="\t" '{ if ($3!="aou") $3="oth"; print $1, $idx"_"$3 }' \
| cat <( echo -e "sample_id\tlabel" ) - \
> gatksv_05B_outliers/dfci-g2c.intake_pop_labels.aou_split.tsv
for alg in depth manta melt wham; do
  niqr=6
  case $alg in
    "depth")
      niqr=10
      prefix="Depth"
      ;;
    "manta")
      prefix="Manta"
      ;;
    "melt")
      prefix="Melt"
      ;;
    "wham")
      prefix="Wham"
      ;;
  esac
  code/scripts/define_variant_count_outlier_samples.R \
    --counts-tsv gatksv_05B_outliers/dfci-g2c.05_sv_counts.$alg.tsv \
    --sample-labels-tsv gatksv_05B_outliers/dfci-g2c.intake_pop_labels.aou_split.tsv \
    --n-iqr $niqr \
    --no-lower-filter \
    --plot \
    --plot-title-prefix $prefix \
    --out-prefix gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.$alg
done

# Merge outlier sample lists
cat gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.*.outliers.samples.list \
| sort -V | uniq \
> gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.all_outlier_samples.list

# Update sample metadata with 05B outlier failure labels
code/scripts/append_qc_fail_metadata.R \
  --qc-tsv dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz \
  --new-column-name clusterbatch_qc_pass \
  --all-samples-list <( cat batch_info/sample_lists/* ) \
  --fail-samples-list gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.all_outlier_samples.list \
  --outfile dfci-g2c.sample_meta.post_clusterbatch.tsv
gzip -f dfci-g2c.sample_meta.post_clusterbatch.tsv

# Compress and archive outlier data for future reference
tar -czvf gatksv_05B_outliers.tar.gz gatksv_05B_outliers
gsutil -m cp \
  gatksv_05B_outliers.tar.gz \
  gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.all_outlier_samples.list \
  dfci-g2c.sample_meta.post_clusterbatch.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/

# Replot sample QC after excluding outliers above
qcplotdir=dfci-g2c.phase1.clusterbatch_qc_pass.plots
if [ ! -e $qcplotdir ]; then mkdir $qcplotdir; fi
code/scripts/plot_intake_qc.R \
  --qc-tsv dfci-g2c.sample_meta.post_clusterbatch.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --pass-column clusterbatch_qc_pass \
  --out-prefix $qcplotdir/dfci-g2c.phase1.clusterbatch_qc_pass
tar -czvf dfci-g2c.phase1.clusterbatch_qc_pass.plots.tar.gz $qcplotdir
gsutil -m cp \
  dfci-g2c.phase1.clusterbatch_qc_pass.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/gatksv_qc/

# Launch workflows to remove outlier samples from clustered VCFs for each batch
module_submission_routine_all_batches 05B


########################
# 05C | ReclusterBatch #
########################

# Note: this is not a canonical GATK-SV module and was instituted specifically for G2C

module_submission_routine_all_batches 05C


#############################
# 06 | GenerateBatchMetrics #
#############################

module_submission_routine_all_batches 06


#########################
# 07 | FilterBatchSites #
#########################

module_submission_routine_all_batches 07


###########################
# 08 | FilterBatchSamples #
###########################

# Note: similar to 05B above, this step uses custom outlier definitions
# The below code needs to be run just once in a single workspace

# This also requires all packages installed and outputs as for 05B above
# If necessary, rerun the 05B library installation steps before proceeding

# Collect SV count data for each sample per algorithm
if ! [ -e data/gatksv_08_outliers ]; then
  mkdir gatksv_08_outliers
fi
for alg in depth manta melt wham; do
  gsutil -m cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/07/*/PlotSVCountsPerSample/CountSVsPerSamplePerType/**.$alg.with_evidence.svcounts.txt \
  | sort -Vrk1,1 -k2,2V | uniq \
  > gatksv_08_outliers/dfci-g2c.08_sv_counts.$alg.tsv
done

# Define outlier samples cohort-wide as Q3 + 6 IQR for all algorithms
for alg in depth manta melt wham; do
  case $alg in
    "depth")
      prefix="Depth"
      ;;
    "manta")
      prefix="Manta"
      ;;
    "melt")
      prefix="Melt"
      ;;
    "wham")
      prefix="Wham"
      ;;
  esac
  code/scripts/define_variant_count_outlier_samples.R \
    --counts-tsv gatksv_08_outliers/dfci-g2c.08_sv_counts.$alg.tsv \
    --sample-labels-tsv gatksv_05B_outliers/dfci-g2c.intake_pop_labels.aou_split.tsv \
    --n-iqr 6 \
    --no-lower-filter \
    --plot \
    --plot-title-prefix $prefix \
    --out-prefix gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.$alg
done

# Merge outlier sample lists
cat gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.*.outliers.samples.list \
| sort -V | uniq \
> gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.all_outlier_samples.list

# Update sample metadata with 08 outlier failure labels
code/scripts/append_qc_fail_metadata.R \
  --qc-tsv dfci-g2c.sample_meta.post_clusterbatch.tsv.gz \
  --new-column-name filtersites_qc_pass \
  --all-samples-list <( cat batch_info/sample_lists/* | fgrep -wvf \
                        gatksv_05B_outliers/dfci-g2c.gatksv.05B_outliers.all_outlier_samples.list ) \
  --fail-samples-list gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.all_outlier_samples.list \
  --outfile dfci-g2c.sample_meta.post_filtersites.tsv
gzip -f dfci-g2c.sample_meta.post_filtersites.tsv

# Compress and archive outlier data for future reference
tar -czvf gatksv_08_outliers.tar.gz gatksv_08_outliers
gsutil -m cp \
  gatksv_08_outliers.tar.gz \
  gatksv_08_outliers/dfci-g2c.gatksv.08_outliers.all_outlier_samples.list \
  dfci-g2c.sample_meta.post_filtersites.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/

# Replot sample QC after excluding outliers above
qcplotdir=dfci-g2c.phase1.filtersites_qc_pass.plots
if [ ! -e $qcplotdir ]; then mkdir $qcplotdir; fi
code/scripts/plot_intake_qc.R \
  --qc-tsv dfci-g2c.sample_meta.post_filtersites.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --pass-column clusterbatch_qc_pass \
  --pass-column filtersites_qc_pass \
  --out-prefix $qcplotdir/dfci-g2c.phase1.filtersites_qc_pass
tar -czvf dfci-g2c.phase1.filtersites_qc_pass.plots.tar.gz $qcplotdir
gsutil -m cp \
  dfci-g2c.phase1.filtersites_qc_pass.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/gatksv_qc/

# Launch workflows to remove outlier samples from filtered VCFs for each batch
module_submission_routine_all_batches 08


########################
# 09 | MergeBatchSites #
########################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Submit workflow
submit_cohort_module 09

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.09-MergeBatchSites.job_ids.list )

# Once complete, stage outputs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.09-MergeBatchSites.job_ids.list ) \
| awk '{ print $2 }' | gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/09/

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell/*/MergeBatchSites/** >> uris_to_delete.list
cleanup_garbage


######################
# 10 | GenotypeBatch #
######################

module_submission_routine_all_batches 10


#######################
# 11 | RegenotypeCNVs #
#######################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Submit workflow
submit_cohort_module 11

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.11-RegenotypeCNVs.job_ids.list )

# Once complete, stage outputs and cleanup garbage
stage_cohort_module 11
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.11-RegenotypeCNVs.job_ids.list ) \
| awk '{ print $2 }' | gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/11/

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell/*/RegenotypeCNVs/** >> uris_to_delete.list
cleanup_garbage


#######################
# 12 | CombineBatches #
#######################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Note 2: this module is handled differently by submit_cohort_module since it's
# parallelized by chromosome with 24 independent submissions

# All cleanup and tracking is handled by a helper routine within submit_cohort_module

submit_cohort_module 12


###############################
# 13 | ResolveComplexVariants #
###############################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Submit workflow
submit_cohort_module 13

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.13-ResolveComplexVariants.job_ids.list )

# # Once complete, stage outputs and cleanup garbage
# stage_cohort_module 11
# cromshell -t 120 --no_turtle -mc list-outputs \
#   $( tail -n1 cromshell/job_ids/dfci-g2c.v1.11-RegenotypeCNVs.job_ids.list ) \
# | awk '{ print $2 }' | gsutil -m cp -I \
#   $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/11/

# # Once staged, clean up outputs
# gsutil -m ls $WORKSPACE_BUCKET/cromwell/*/RegenotypeCNVs/** >> uris_to_delete.list
# cleanup_garbage
