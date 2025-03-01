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


#####################
# Prepare intervals #
#####################

# Note: this only needs to be run once across all workspaces

# Refresh staging directory
staging_dir=staging/PrepIntervals
if [ -e $staging_dir ]; then rm -rf $staging_dir; fi; mkdir $staging_dir

# Download Broad standard hg38 calling intervals list
gsutil -m cp \
  gs://gcp-public-data--broad-references/hg38/v0/wgs_calling_regions.hg38.interval_list \
  $staging_dir/

# Split intervals into 24 primary chromosomes
ilist=$staging_dir/wgs_calling_regions.hg38.interval_list
for k in $( seq 1 22 ) X Y; do
  head -n1 $ilist > $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.interval_list
  fgrep -w "SN:chr$k" $ilist \
  >> $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.interval_list
  awk -v contig="chr$k" '{ if ($1==contig) print }' $ilist \
  >> $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.interval_list
done

# Estimate total genomic space to be processed by each of five workspaces
cat $staging_dir/gatkhc.wgs_calling_regions.hg38.chr*.interval_list \
| fgrep -v "@" | cut -f1-3 | bedtools merge -i - \
| awk '{ sum+=$3-$2 }END{ printf "%i\n", sum / 5 }'

# Subdivide intervals per chromosome to allow maximal parallelization
# Logic is as follows: ~600Mb per workspace, 
# Logic as follows: ~600Mb per workspace / 2900 task quota = 200kb avg interval size
for k in $( seq 1 22 ) X Y; do
  code/scripts/split_intervals.py \
    -i $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.interval_list \
    -t 200000 \
    --gatk-style \
    -o $staging_dir/gatkhc.wgs_calling_regions.hg38.chr$k.sharded.intervals
done

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
for k in $( seq 1 22 ) X Y; do
  kc=$( fgrep -v "@" \
          staging/PrepIntervals/gatkhc.wgs_calling_regions.hg38.chr$k.sharded.intervals \
        | wc -l | awk '{ printf "%i\n", $1 }' )
  echo "\"chr$k\" : {\"CONTIG_SCATTER_COUNT\" : $kc },"
done | paste -s -d\  | sed 's/,$//g' \
>> $staging_dir/contig_variable_overrides.json
echo " }" >> $staging_dir/contig_variable_overrides.json

# Write template .json for input
cat << EOF > $staging_dir/GnarlyJointGenotypingPart1.inputs.template.json
{
  "GnarlyJointGenotypingPart1.callset_name": "dfci-g2c.v1.\$CONTIG",
  "GnarlyJointGenotypingPart1.dbsnp_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
  "GnarlyJointGenotypingPart1.GnarlyGenotyper.disk_size_gb": 100,
  "GnarlyJointGenotypingPart1.GnarlyGenotyper.machine_mem_mb": 20000,
  "GnarlyJointGenotypingPart1.gnarly_scatter_count": 1,
  "GnarlyJointGenotypingPart1.import_gvcf_batch_size": 100,
  "GnarlyJointGenotypingPart1.ImportGVCFs.machine_mem_mb": 48000,
  "GnarlyJointGenotypingPart1.make_hard_filtered_sites": false,
  "GnarlyJointGenotypingPart1.medium_disk": 125,
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
  --max-attempts 4




#### DEV - carefully benchmarking & optimizing resources just on chr19

export CONTIG=chr19
export CONTIG_SCATTER_COUNT=279

code/scripts/envsubst.py \
  -i $staging_dir/GnarlyJointGenotypingPart1.inputs.template.json \
  -o cromshell/inputs/GnarlyJointGenotypingPart1.inputs.chr19.json

cromshell --no_turtle -t 120 -mc submit \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  code/wdl/gatk-hc/GnarlyJointGenotypingPart1.wdl \
  cromshell/inputs/GnarlyJointGenotypingPart1.inputs.chr19.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.JointGenotyping.chr19.job_ids.list

while [ TRUE ]; do
  echo -e "\n\n"
  date
  cromshell -t 120 counts -x \
    $( tail -n1 cromshell/job_ids/dfci-g2c.v1.JointGenotyping.chr19.job_ids.list )
  sleep 10m
done

# Optimized resources based on chr19 -- need to update input json above once workflow completes
# import_gvcf_disk_gb = 40
# "GnarlyJointGenotypingPart1.GnarlyGenotyper.disk_size_gb": 50,
# "GnarlyJointGenotypingPart1.GnarlyGenotyper.machine_mem_mb": 20,



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

