#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-HC joint genotyping on G2C phase 1

# Note that this code is designed to be run inside the AoU Researcher Workbench


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in staging cromshell cromshell/inputs cromshell/job_ids cromshell/progress; do
  if ! [ -e $dir ]; then
    mkdir $dir
  fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
find code/ -name "*.py" | xargs -I {} chmod a+x {}
find code/ -name "*.R" | xargs -I {} chmod a+x {}

# Source .bashrc and bash utility functions
. code/refs/dotfiles/aou.rw.bashrc
. code/refs/general_bash_utils.sh

# Install necessary packages
. code/refs/install_packages.sh python

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Create dependencies .zip for workflow submissions
cd code/wdl/pancan_germline_wgs && \
zip g2c.dependencies.zip *.wdl && \
mv g2c.dependencies.zip ~/ && \
cd ~

# Infer workspace number and save as environment variable
export WN=$( get_workspace_number )

# Download workspace-specific contig lists
gsutil cp -r \
  gs://dfci-g2c-refs/hg38/contig_lists \
  ./


#######################
# Generate sample map #
#######################

# Note: this only needs to be run once across all workspaces

# Refresh staging directory
staging_dir=staging/GenerateSampleMap
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Gather information needed for all samples to be included for joint genotyping
gsutil -m cat \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/qc-filtering/dfci-g2c.sample_meta.posthoc_outliers.tsv.gz \
| gunzip -c \
| awk -v FS="\t" -v OFS="\t" '{ if ($NF=="True") print $1, $2, $3 }' \
| sort -Vk1,1 \
> $staging_dir/dfci-g2c.v1.gatkhc.input_samples.info.tsv

# Write GATK-HC sample map
while read sid oid cohort; do
  if [ $cohort == "aou" ]; then
    uri_prefix=$MAIN_WORKSPACE_BUCKET
  else
    uri_prefix="gs:/"
  fi
  echo -e "$sid\t$uri_prefix/dfci-g2c-inputs/$cohort/gatk-hc/reblocked/$oid.reblocked.g.vcf.gz"
done < $staging_dir/dfci-g2c.v1.gatkhc.input_samples.info.tsv \
> $staging_dir/dfci-g2c.v1.gatkhc.sample_map.tsv

# Move sample map to main bucket for Cromwell access
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.gatkhc.sample_map.tsv \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/


########################################################################
# Collect summary distributions from gnomAD v4.1 for interval sharding #
########################################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Refresh staging directory
staging_dir=staging/GnomadSiteMetrics
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Write template .json of inputs for chromsharded manager
cat << EOF > $staging_dir/PreprocessGnomadSiteMetrics.inputs.template.json
{
  "PreprocessGnomadSiteMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:initial_mar3",
  "PreprocessGnomadSiteMetrics.g2c_pipeline_docker": "vanallenlab/g2c_pipeline:sv_counting",
  "PreprocessGnomadSiteMetrics.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "PreprocessGnomadSiteMetrics.min_af_bin": 0.00003618403,
  "PreprocessGnomadSiteMetrics.output_prefix": "gnomad.v4.1.\$CONTIG",
  "PreprocessGnomadSiteMetrics.snv_n_samples": 76215,
  "PreprocessGnomadSiteMetrics.snv_scatter_intervals": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/gatkhc.wgs_calling_regions.hg38.\$CONTIG.sharded.intervals",
  "PreprocessGnomadSiteMetrics.snv_vcf": "gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.\$CONTIG.vcf.bgz",
  "PreprocessGnomadSiteMetrics.snv_vcf_idx": "gs://gcp-public-data--gnomad/release/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.\$CONTIG.vcf.bgz.tbi",
  "PreprocessGnomadSiteMetrics.sv_n_samples": 63046,
  "PreprocessGnomadSiteMetrics.sv_scatter_intervals": "gs://dfci-g2c-refs/hg38/contig_lists/\$CONTIG.list",
  "PreprocessGnomadSiteMetrics.sv_vcf": "gs://gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz",
  "PreprocessGnomadSiteMetrics.sv_vcf_idx": "gs://gcp-public-data--gnomad/release/4.1/genome_sv/gnomad.v4.1.sv.sites.vcf.gz.tbi"
}
EOF

# Submit, monitor, stage, and cleanup short variant QC metadata task
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/PreprocessGnomadSiteMetrics.wdl \
  --input-json-template $staging_dir/PreprocessGnomadSiteMetrics.inputs.template.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/refs/gnomad_v4_site_metrics \
  --name GnomadSiteMetrics \
  --status-tsv cromshell/progress/dfci-g2c.v1.initial_qc.GnomadSiteMetrics.progress.tsv \
  --workflow-id-log-prefix "gnomad.v4.site_metrics" \
  --outer-gate 30 \
  --submission-gate 30 \
  --max-attempts 2


#####################
# Prepare intervals #
#####################

# Refresh staging directory
staging_dir=staging/PrepIntervals
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Download Broad standard hg38 calling intervals list
gsutil -m cp \
  gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list \
  $staging_dir/

# Split intervals into one per primary chromosome
ilist=$staging_dir/wgs_calling_regions.hg38.interval_list
for k in $( seq 1 22 ) X Y; do
  head -n1 $ilist > $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.interval_list
  fgrep -w "SN:chr$k" $ilist \
  >> $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.interval_list
  awk -v contig="chr$k" '{ if ($1==contig) print }' $ilist \
  >> $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.interval_list
done

# Estimate ideal number of shards per chromosome
# Math as follows: 2,900 task quota x 5 workspaces = 14,500 total shards
# We will increase this value by 50% to ensure we aren't wasting much quota
# without also overloading the cromwell server
# Shards per chrom = 14,500 x 1.5 x chromosome size / genome size
genome_size=$( cat $staging_dir/gatkhc.wgs_calling_regions.hg38.chr*.interval_list \
               | fgrep -v "@" | cut -f1-3 | bedtools merge -i - \
               | awk '{ sum+=$3-$2 }END{ printf "%i\n", sum }' )

# Re-shard GATK-HC intervals per chromosome
while read contig; do
  # Get size of contig
  contig_size=$( cat $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.interval_list \
                 | fgrep -v "@" | cut -f1-3 | bedtools merge -i - \
                 | awk '{ sum+=$3-$2 }END{ printf "%i\n", sum }' )

  # Calculate number of desired shards
  n_shards=$( echo "14500" \
              | awk -v denom=$genome_size -v numer=$contig_size \
                '{ printf "%i\n", $1 * 1.5 * numer / denom }' )

  # Download gnomAD variant sites from this contig
  gsutil -m cp \
    $MAIN_WORKSPACE_BUCKET/refs/gnomad_v4_site_metrics/$contig/*/gnomad.v4.1.$contig.*.sites.bed.gz* \
    $staging_dir/
  find $staging_dir -name "*bed.gz" | xargs -I {} tabix -f {}

  # Estimate ideal number of variants per shard
  n_vars=$( zcat $staging_dir/gnomad.v4.1.$contig.*.sites.bed.gz | wc -l | awk '{ print $1-3 }' )
  vars_per_shard=$( echo $n_vars | awk -v denom=$n_shards '{ printf "%i\n", $1 / denom }' )

  # Shard GATK-HC intervals according to parameters determined above
  code/scripts/split_intervals.py \
    -i $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.interval_list \
    --var-sites $staging_dir/gnomad.v4.1.$contig.snv.sites.bed.gz \
    --var-sites $staging_dir/gnomad.v4.1.$contig.indel.sites.bed.gz \
    --var-sites $staging_dir/gnomad.v4.1.$contig.sv.sites.bed.gz \
    --vars-per-shard $vars_per_shard \
    --gatk-style \
    --verbose \
    -o $staging_dir/gatkhc.wgs_calling_regions.hg38.$contig.sharded.intervals

  # Clean up
  rm $staging_dir/gnomad.v4.1.*.*.sites.bed.gz*

done < contig_lists/dfci-g2c.v1.contigs.w$WN.list

# Copy all sharded interval lists to google bucket for Cromwell access
gsutil cp \
  $staging_dir/gatkhc.wgs_calling_regions.hg38.chr*.sharded.intervals \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/


####################
# Joint genotyping #
####################

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Refresh staging directory
staging_dir=staging/JointGenotyping
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Check to ensure there is a local copy of calling intervals
if ! [ -e staging/PrepIntervals ]; then
  mkdir staging/PrepIntervals
  gsutil -m cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/*.sharded.intervals \
    staging/PrepIntervals/
fi

# Write .json of contig-specific scatter counts
echo "{ " > $staging_dir/contig_variable_overrides.json
while read contig; do
  kc=$( fgrep -v "@" \
          staging/PrepIntervals/gatkhc.wgs_calling_regions.hg38.$contig.sharded.intervals \
        | wc -l | awk '{ printf "%i\n", $1 }' )
  echo "\"$contig\" : {\"CONTIG_SCATTER_COUNT\" : $kc },"
done < contig_lists/dfci-g2c.v1.contigs.w$WN.list \
| paste -s -d\  | sed 's/,$//g' \
>> $staging_dir/contig_variable_overrides.json
echo " }" >> $staging_dir/contig_variable_overrides.json

# Write template .json for input
cat << EOF > $staging_dir/GnarlyJointGenotypingPart1.inputs.template.json
{
  "GnarlyJointGenotypingPart1.callset_name": "dfci-g2c.v1.\$CONTIG",
  "GnarlyJointGenotypingPart1.dbsnp_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
  "GnarlyJointGenotypingPart1.GnarlyGenotyper.disk_size_gb": 50,
  "GnarlyJointGenotypingPart1.GnarlyGenotyper.machine_mem_mb": 20000,
  "GnarlyJointGenotypingPart1.gnarly_scatter_count": 1,
  "GnarlyJointGenotypingPart1.import_gvcfs_batch_size": 100,
  "GnarlyJointGenotypingPart1.import_gvcfs_disk_gb": 40,
  "GnarlyJointGenotypingPart1.ImportGVCFs.machine_mem_mb": 48000,
  "GnarlyJointGenotypingPart1.intervals_already_split": true,
  "GnarlyJointGenotypingPart1.make_hard_filtered_sites": false,
  "GnarlyJointGenotypingPart1.ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
  "GnarlyJointGenotypingPart1.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "GnarlyJointGenotypingPart1.ref_fasta_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai",
  "GnarlyJointGenotypingPart1.sample_name_map": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/dfci-g2c.v1.gatkhc.sample_map.tsv",
  "GnarlyJointGenotypingPart1.top_level_scatter_count": \$CONTIG_SCATTER_COUNT,
  "GnarlyJointGenotypingPart1.unpadded_intervals_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/gatkhc.wgs_calling_regions.hg38.\$CONTIG.sharded.intervals"
}
EOF

# Joint genotype per chromosome using chromsharded manager
# Note: all cleanup and tracking is handled by the chromshard manager, so
# no outputs or execution buckets need to be manually staged/cleared
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-hc/GnarlyJointGenotypingPart1.wdl \
  --input-json-template $staging_dir/GnarlyJointGenotypingPart1.inputs.template.json \
  --contig-variable-overrides $staging_dir/contig_variable_overrides.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/ \
  --name JointGenotyping \
  --contig-list contig_lists/dfci-g2c.v1.contigs.w$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.JointGenotyping.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 60 \
  --submission-gate 5 \
  --max-attempts 2


###############
# VCF cleanup #
###############

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Refresh staging directory
staging_dir=staging/PosthocCleanup
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Write template .json for input
cat << EOF > $staging_dir/PosthocCleanupPart1.inputs.template.json
{
  "PosthocCleanupPart1.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PosthocCleanupPart1.g2c_pipeline_docker": "vanallenlab/g2c_pipeline:sv_counting",
  "PosthocCleanupPart1.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "PosthocCleanupPart1.ref_fasta": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta",
  "PosthocCleanupPart1.unpadded_intervals_file": "$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/gatkhc.wgs_calling_regions.hg38.\$CONTIG.sharded.intervals",
  "PosthocCleanupPart1.vcf": "TBD",
  "PosthocCleanupPart1.vcf_idx": "TBD"
}
EOF

# Perform post hoc VCF cleanup (split multiallelics, minimize indel representation)
# Also count variants by type per sample (needed for outlier definition; see below)
code/scripts/manage_chromshards.py \
  --wdl code/wdl/gatk-hc/PosthocCleanupPart1.wdl \
  --input-json-template $staging_dir/PosthocCleanupPart1.inputs.template.json \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart1/ \
  --contig-list contig_lists/dfci-g2c.v1.contigs.w$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.PosthocCleanupPart1.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 45 \
  --max-attempts 3


###########################
# Outlier sample analysis #
###########################

# Note that this section only needs to be run from one workspace for the entire cohort

# Reaffirm staging directory
staging_dir=staging/PosthocCleanup

# TODO: implement this


#############################################################
# Exclude outlier samples and apply site-level hard filters #
#############################################################

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Reaffirm staging directory
staging_dir=staging/PosthocCleanup

# TODO: implement this


###########################################
# Exclude outliers also from GATK-SV VCFs #
###########################################

# Note that this section only needs to be run from one workspace for the entire cohort

# Reaffirm staging directory
staging_dir=staging/PosthocCleanup

# TODO: implement this

