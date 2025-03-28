#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Quality control and filtering of G2C germline callset after joint genotyping

# Note that this code is designed to be run inside the AoU Researcher Workbench


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in cromshell cromshell/inputs cromshell/inputs/templates \
           cromshell/job_ids cromshell/progress staging; do
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
. code/refs/install_packages.sh python R

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


###########################################
# Compile .fam file of reported relatives #
###########################################

# Initialize staging directory
staging_dir=staging/fam_curation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Download sample metadata after GATK-HC
gsutil cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/

# Localize external .fam files staged in G2C bucket
gsutil cp -r \
  gs://dfci-g2c-refs/fam \
  $staging_dir/

# Curate 1000 Genomes project trios
code/scripts/external_fam_to_g2c.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --in-fam $staging_dir/fam/hgsvc.ped \
  --cohort hgsvc \
  --out-fam $staging_dir/dfci-g2c.hgsvc.fam

# Curate CEPH pedigrees
code/scripts/external_fam_to_g2c.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --in-fam $staging_dir/fam/ceph.remapped_ids.fam \
  --cohort ceph \
  --out-fam $staging_dir/dfci-g2c.ceph.fam

# Curate UFC pedigrees
code/scripts/external_fam_to_g2c.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --in-fam $staging_dir/fam/dfci_ufc_cases.fam \
  --cohort ufc \
  --out-fam $staging_dir/dfci-g2c.ufc.fam

# Curate MESA pedigrees
code/scripts/external_fam_to_g2c.R \
  --metadata-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --in-fam $staging_dir/fam/mesa.fam \
  --cohort mesa \
  --out-fam $staging_dir/dfci-g2c.mesa.fam

# Collapse all cohort-specific .fam files
cat $staging_dir/dfci-g2c.*.fam \
| sort -Vk1,1 -k2,2V \
> $staging_dir/dfci-g2c.reported_families.fam

# Stage combined .fam file in Google bucket for future use
gsutil cp \
  $staging_dir/dfci-g2c.reported_families.fam \
  $MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/


############################################
# Collect initial short variant QC metrics #
############################################

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Build chromosome-specific override json of VCFs and VCF indexes
add_contig_vcfs_to_chromshard_overrides_json \
  $staging_dir/CollectShortVariantQcMetrics.contig_variable_overrides.json \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart2 \
  filtered_vcfs \
  filtered_vcf_idxs

# Write template input .json for short variant QC metric collection
cat << EOF > $staging_dir/CollectShortVariantQcMetrics.inputs.template.json
{
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.common_af_cutoff": 0.001,
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:f046423",
  "CollectVcfQcMetrics.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 1000,
  "CollectVcfQcMetrics.output_prefix": "dfci-g2c.v1.gatkhc.initial_qc.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 7.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.trios_fam_file": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.vcfs": \$CONTIG_VCFS,
  "CollectVcfQcMetrics.vcf_idxs": \$CONTIG_VCF_IDXS
}
EOF

# Submit, monitor, stage, and cleanup short variant QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectShortVariantQcMetrics.inputs.template.json \
  --contig-variable-overrides $staging_dir/CollectShortVariantQcMetrics.contig_variable_overrides.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/ShortVariantMetrics/ \
  --name CollectInitialShortVariantQcMetrics \
  --contig-list contig_lists/dfci-g2c.v1.contigs.w$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.CollectShortVariantQcMetrics.initial_qc.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --submission-gate 10 \
  --max-attempts 2


#################################
# Collect initial SV QC metrics #
#################################

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Write template input .json for short variant QC metric collection
cat << EOF > $staging_dir/CollectSVQcMetrics.inputs.template.json
{
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.common_af_cutoff": 0.001,
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:26ce17a",
  "CollectVcfQcMetrics.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 1000,
  "CollectVcfQcMetrics.n_records_per_shard": 10000,
  "CollectVcfQcMetrics.output_prefix": "dfci-g2c.v1.gatksv.initial_qc.\$CONTIG",
  "CollectVcfQcMetrics.shard_vcf": true,
  "CollectVcfQcMetrics.trios_fam_file": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.vcfs": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset/\$CONTIG/HardFilterPart2/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.posthoc_filtered.vcf.gz"],
  "CollectVcfQcMetrics.vcf_idxs": ["$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset/\$CONTIG/HardFilterPart2/dfci-g2c.v1.\$CONTIG.concordance.gq_recalibrated.posthoc_filtered.vcf.gz.tbi"]
}
EOF

# Submit, monitor, stage, and cleanup SV QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectSVQcMetrics.inputs.template.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/SVMetrics/ \
  --name CollectInitialSvQcMetrics \
  --contig-list contig_lists/dfci-g2c.v1.contigs.w$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.CollectSvQcMetrics.initial_qc.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --outer-gate 30 \
  --max-attempts 2


##########################################
# Analyze & visualize initial QC metrics #
##########################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Build input arrays
for key in size_distrib af_distrib size_vs_af_distrib \
              all_svs_bed common_snvs_bed common_indels_bed common_svs_bed; do
  if [ -e $staging_dir/$key.uris.list ]; then
    rm $staging_dir/$key.uris.list
  fi
done
for k in $( seq 1 22 ) X Y; do
  
  # Localize output tracker jsons
  gsutil cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/ShortVariantMetrics/chr$k/CollectShortVariantQcMetrics.chr$k.outputs.json \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/SVMetrics/chr$k/CollectSvQcMetrics.chr$k.outputs.json \
    $staging_dir/

  # Get URIs for QC metrics
  for jprefix in CollectShortVariantQcMetrics CollectSvQcMetrics; do
    json=$staging_dir/$jprefix.chr$k.outputs.json
    for key in size_distrib af_distrib size_vs_af_distrib \
               all_svs_bed common_snvs_bed common_indels_bed common_svs_bed; do
      jq .$key $json | fgrep -xv "null" | tr -d '"' >> $staging_dir/$key.uris.list
    done
  done

  # Clear local copies of output tracker jsons
  rm \
    $staging_dir/CollectShortVariantQcMetrics.chr$k.outputs.json \
    $staging_dir/CollectSvQcMetrics.chr$k.outputs.json
done

# Write input .json for short variant QC metric collection
cat << EOF > cromshell/inputs/PlotInitialVcfQcMetrics.inputs.json
{
  "PlotVcfQcMetrics.af_distribution_tsvs": $( collapse_txt $staging_dir/af_distrib.uris.list ),
  "PlotVcfQcMetrics.all_sv_beds": $( collapse_txt $staging_dir/all_svs_bed.uris.list ),
  "PlotVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PlotVcfQcMetrics.common_af_cutoff": 0.001,
  "PlotVcfQcMetrics.common_snv_beds": $( collapse_txt $staging_dir/common_snvs_bed.uris.list ),
  "PlotVcfQcMetrics.common_indel_beds": $( collapse_txt $staging_dir/common_indels_bed.uris.list ),
  "PlotVcfQcMetrics.common_sv_beds": $( collapse_txt $staging_dir/common_svs_bed.uris.list ),
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:26ce17a",
  "PlotVcfQcMetrics.output_prefix": "dfci-g2c.v1.gatksv.initial_qc",
  "PlotVcfQcMetrics.size_distribution_tsvs": $( collapse_txt $staging_dir/size_distrib.uris.list ),
  "PlotVcfQcMetrics.size_vs_af_distribution_tsvs": $( collapse_txt $staging_dir/size_vs_af_distrib.uris.list )
}
EOF

# Submit QC visualization workflow
cromshell --no_turtle -t 120 -mc submit \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip g2c.dependencies.zip \
  code/wdl/pancan_germline_wgs/vcf-qc/PlotVcfQcMetrics.wdl \
  cromshell/inputs/PlotInitialVcfQcMetrics.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.PlotInitialVcfQcMetrics.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotInitialVcfQcMetrics.job_ids.list )

# # Once patches are complete, manually stage output 
# patch_wid=$( tail -n1 cromshell/job_ids/GnarlyJointGenotypingPart1.inputs.$contig.patch.job_ids.list )
# gsutil -m cp \
#   $WORKSPACE_BUCKET/cromwell/execution/GnarlyJointGenotypingPart1/$patch_wid/**/call-GnarlyGenotyperFT/**dfci-g2c.v1.$contig.*.vcf.gz* \
#   $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/JointGenotyping/$contig/

# # Clear Cromwell execution & output buckets for patch jobs
# gsutil -m ls $( cat cromshell/job_ids/GnarlyJointGenotypingPart1.inputs.$contig.patch.job_ids.list \
#                 | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell/*/GnarlyJointGenotypingPart1/" \
#                   '{ print bucket_prefix$1"/**" }' ) \
# > uris_to_delete.list
# cleanup_garbage

