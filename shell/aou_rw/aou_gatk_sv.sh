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

# Once complete, stage outputs
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

# Once complete, stage outputs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.13-ResolveComplexVariants.job_ids.list ) \
| awk '{ print $2 }' | gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/13/

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell/*/ResolveComplexVariants/** >> uris_to_delete.list
cleanup_garbage


###############################
# 14A | FilterCoverageSamples #
###############################

# Note: this is not a canonical GATK-SV module and was instituted specifically for G2C
# All this does is remove outliers from the median coverage file to avoid 
# RDtest normalization errors in module 14

# Get list of samples present in VCF at this stage
if [ ! -e staging/14A-FilterCoverageSamples ]; then
  mkdir staging/14A-FilterCoverageSamples
fi
gsutil cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/13/dfci-g2c.v1.reshard_vcf.chrY.resharded.vcf.gz \
| gunzip -c | head -n300 | fgrep -v "##" | fgrep "#" | cut -f10- \
| sed 's/\t/\n/g' | sort -V | uniq \
> staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list
gsutil -m cp \
  staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/

module_submission_routine_all_batches 14A


################################
# 14 | GenotypeComplexVariants #
################################

# Note: this module only needs to be run once in one workspace for the whole cohort

# # All cleanup and tracking is handled by a helper routine within submit_cohort_module

# Subset .ped file to those present in VCF
samples14A_list=staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list
if [ ! -e $samples14A_list ]; then
  echo "ERROR: cannot locate $samples14A_list. Regenerate if necessary."
else
  if [ ! -e staging/14-GenotypeComplexVariants ]; then
    mkdir staging/14-GenotypeComplexVariants
  fi
  gsutil -m cat \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/dfci-g2c.all_samples.ped \
  | fgrep -wf staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list \
  > staging/14-GenotypeComplexVariants/dfci-g2c.all_samples.gatksv_module14.ped
  gsutil cp \
    staging/14-GenotypeComplexVariants/dfci-g2c.all_samples.gatksv_module14.ped \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/
fi

# Submit workflow
submit_cohort_module 14

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.14-GenotypeComplexVariants.job_ids.list )

# Once complete, stage outputs
cromshell -t 120 --no_turtle -mc list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.14-GenotypeComplexVariants.job_ids.list ) \
| awk '{ print $2 }' | gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/14/

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell/*/GenotypeComplexVariants/** >> uris_to_delete.list
cleanup_garbage


#################
# 15 | CleanVcf #
#################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Submit workflow
submit_cohort_module 15

# Monitor submission
monitor_workflow \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.15-CleanVcf.job_ids.list )

# CleanVcf processes all chromosomes in parallel but is designed to concatenate 
# all chromosomes into a single output VCF. However, since all of our downstream
# steps are parallel by chromosome, there is no point in waiting for CleanVcf to
# finish concatenating the massive overall VCF. Thus, we can determine whether 
# CleanVcf is finished by checking if the concatenation task is present in the
# execution bucket, in which case we can kill the job and relocalize all of the
# chromosome-sharded outputs
cvcf_exec_base=$WORKSPACE_BUCKET/cromwell/execution/CleanVcf/$( tail -n1 cromshell/job_ids/dfci-g2c.v1.15-CleanVcf.job_ids.list )
if [ $( gsutil ls $cvcf_exec_base/call-ConcatCleanedVcfs | wc -l ) -gt 0 ]; then
  for k in $( seq 0 23 ); do
    contig_wid=$( basename $( gsutil ls $cvcf_exec_base/call-CleanVcfChromosome/shard-$k/CleanVcfChromosome/ ) )
    cromshell -t 120 --no_turtle -mc list-outputs $contig_wid \
    | awk '{ print $2 }' | gsutil -m cp -I \
      $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/15/
  done
  cromshell -t 120 abort $( tail -n1 cromshell/job_ids/dfci-g2c.v1.15-CleanVcf.job_ids.list )
fi

# Once staged, clean up outputs
gsutil -m ls $WORKSPACE_BUCKET/cromwell/*/CleanVcf/** >> uris_to_delete.list
cleanup_garbage


##############################
# 16 | RefineComplexVariants #
##############################

# Note: this module only needs to be run once in one workspace for the whole cohort

# Subset batch sample membership files to those present in VCF
samples14A_list=staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list
if [ ! -e $samples14A_list ]; then
  echo "ERROR: cannot locate $samples14A_list. Regenerate if necessary."
else
  for dir in 16-RefineComplexVariants \
             16-RefineComplexVariants/module14A_batch_sample_lists; do
    if [ ! -e staging/$dir ]; then mkdir staging/$dir; fi
  done
  while read bid; do
    fgrep -xf $samples14A_list batch_info/sample_lists/$bid.samples.list \
    > staging/16-RefineComplexVariants/module14A_batch_sample_lists/$bid.gatksv_module14A.samples.list
  done < batch_info/dfci-g2c.gatk-sv.batches.list
  gsutil cp -r \
    staging/16-RefineComplexVariants/module14A_batch_sample_lists \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/refs/
fi

# Note 2: this module is handled differently by submit_cohort_module since it's
# parallelized by chromosome with 24 independent submissions

# All cleanup and tracking is handled by a helper routine within submit_cohort_module

submit_cohort_module 16


##############################################################
# Post hoc outlier sample exclusion and hard variant filters #
##############################################################

# Write input .json for SV counting task
staging_dir=staging/posthoc_filter
if [ ! -e $staging_dir ]; then mkdir $staging_dir; fi
for k in $( seq 1 22 ) X Y; do
  # TODO: write path to all 24 chromsharded VCFs to temp file 
  # >> $staging_dir/vcfs.list
  # TODO: write path to all 24 chromsharded VCF indexes to temp file 
  # >> $staging_dir/vcf_idxs.list
done
cat << EOF > cromshell/inputs/count_svs_posthoc.inputs.json
{
  "CountSvsPerSample.g2c_pipeline_docker": "String",
  "CountSvsPerSample.output_prefix": "dfci-g2c.v1.gatksv_postCleanVcf",
  "CountSvsPerSample.sv_pipeline_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-pipeline:2025-01-14-v1.0.1-88dbd052",
  "CountSvsPerSample.vcfs": $( collapse_txt $staging_dir/vcfs.list ),
  "CountSvsPerSample.vcf_idxs": $( collapse_txt $staging_dir/vcf_idxs.list )
}
EOF

# Submit SV counting task
cromshell --no_turtle -t 120 -mc submit \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  code/wdl/gatk-sv/CountSvsPerSample.wdl \
  cromshell/inputs/count_svs_posthoc.inputs.json
| jq .id | tr -d '"' \
>> cromshell/job_ids/count_svs_posthoc.job_ids.list

# Monitor SV counting task
monitor_workflow \
  $( tail -n1 cromshell/job_ids/count_svs_posthoc.job_ids.list )

# Once complete, download SV counts per sample and exclude CTXs, CNVs, and BNDs
# TODO: implement this

# Define outliers as in 05B or 08 above (maybe don't consider CNV or BND types?)
code/scripts/define_variant_count_outlier_samples.R \
  --counts-tsv gatksv_08_outliers/dfci-g2c.08_sv_counts.$alg.tsv \
  --sample-labels-tsv gatksv_05B_outliers/dfci-g2c.intake_pop_labels.aou_split.tsv \
  --n-iqr 6 \
  --no-lower-filter \
  --plot \
  --plot-title-prefix "GATK-SV" \
  --out-prefix $staging_dir/dfci-g2c.v1.gatksv.posthoc_outliers

# Update sample metadata with 08 outlier failure labels
code/scripts/append_qc_fail_metadata.R \
  --qc-tsv dfci-g2c.sample_meta.post_filtersites.tsv.gz \
  --new-column-name gatksv_posthoc_qc_pass \
  --all-samples-list <( cat batch_info/sample_lists/* | fgrep -wvf \
                        staging/14A-FilterCoverageSamples/dfci-g2c.gatksv.present_at_module14.samples.list ) \
  --fail-samples-list $staging_dir/dfci-g2c.v1.gatksv.posthoc_outliers.samples.list \
  --outfile dfci-g2c.sample_meta.posthoc_outliers.tsv
gzip -f dfci-g2c.sample_meta.posthoc_outliers.tsv

# Compress and archive outlier data for future reference
cd $staging_dir && \
tar -czvf dfci-g2c.v1.gatksv.posthoc_outliers.tar.gz dfci-g2c.v1.gatksv.posthoc_outliers && \
gsutil -m cp \
  dfci-g2c.v1.gatksv.posthoc_outliers.tar.gz \
  $staging_dir/dfci-g2c.v1.gatksv.posthoc_outliers.samples.list \
  dfci-g2c.sample_meta.posthoc_outliers.tsv.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/

# Replot sample QC after excluding outliers above
qcplotdir=dfci-g2c.phase1.gatksv_posthoc_qc_pass.plots
if [ ! -e $qcplotdir ]; then mkdir $qcplotdir; fi
code/scripts/plot_intake_qc.R \
  --qc-tsv dfci-g2c.sample_meta.posthoc_outliers.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --pass-column clusterbatch_qc_pass \
  --pass-column filtersites_qc_pass \
  --pass-column gatksv_posthoc_qc_pass \
  --out-prefix $qcplotdir/dfci-g2c.phase1.gatksv_posthoc_qc_pass
tar -czvf dfci-g2c.phase1.gatksv_posthoc_qc_pass.plots.tar.gz $qcplotdir
gsutil -m cp \
  dfci-g2c.phase1.gatksv_posthoc_qc_pass.plots.tar.gz \
  $MAIN_WORKSPACE_BUCKET/results/gatksv_qc/

# Write template input .json for outlier exclusion & hard filter task
cat << EOF > $staging_dir/PosthocHardFilter.inputs.template.json
{
  "PosthocHardFilter.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PosthocHardFilter.exclude_samples_list": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/dfci-g2c.v1.gatksv.posthoc_outliers.samples.list",
  "PosthocHardFilter.vcf": "TBD",
  "PosthocHardFilter.vcf_idx": "TBD"
}
EOF

# Submit outlier exclusion & hard filter task using chromsharded manager
# Reminder that this manager script handles staging & cleanup too
code/scripts/manage_chromshards.py \
  --wdl gatk-sv/PosthocHardFilter.wdl \
  --input-json-template $staging_dir/PosthocHardFilter.inputs.template.json \
  --staging-bucket $staging_prefix/PosthocHardFilter \
  --name PosthocHardFilter \
  --status-tsv cromshell/progress/dfci-g2c.v1.PosthocHardFilter.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --gate 30 \
  --max-attempts 2


##################
# Raw callset QC #
##################

# TODO: implement this

