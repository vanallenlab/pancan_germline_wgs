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
cp vcf-qc/*.wdl ./ && \
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


##############################
# Collect initial QC metrics #
##############################

# Note: this workflow is scattered across all five workspaces for max parallelization
# It must be submitted as below in each workspace

# Reaffirm staging directory
staging_dir=staging/initial_qc
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Check to ensure there is a local copy of calling intervals
if ! [ -e $staging_dir/calling_intervals ]; then
  mkdir $staging_dir/calling_intervals
  gsutil -m cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/refs/*.sharded.interval_list \
    $staging_dir/calling_intervals/
fi

# Initialize .json of contig-specific overrieds for SV VCF paths and scatter counts
echo "{ " > $staging_dir/CollectVcfQcMetrics.contig_variable_overrides.json
while read contig; do
  kc=$( fgrep -v "@" \
          $staging_dir/calling_intervals/gatkhc.wgs_calling_regions.hg38.$contig.sharded.interval_list \
        | wc -l | awk '{ printf "%i\n", $1 }' )
  echo "\"$contig\" : {\"CONTIG_SCATTER_COUNT\" : $kc,"
  echo "\"CONTIG_VCFS\" : [\"$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset/$contig/HardFilterPart2/dfci-g2c.v1.$contig.concordance.gq_recalibrated.posthoc_filtered.vcf.gz\"],"
  echo "\"CONTIG_VCF_IDXS\" : [\"$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/ExcludeSnvOutliersFromSvCallset/$contig/HardFilterPart2/dfci-g2c.v1.$contig.concordance.gq_recalibrated.posthoc_filtered.vcf.gz.tbi\"] },"
done < contig_lists/dfci-g2c.v1.contigs.w$WN.list \
| paste -s -d\  | sed 's/,$//g' \
>> $staging_dir/CollectVcfQcMetrics.contig_variable_overrides.json
echo " }" >> $staging_dir/CollectVcfQcMetrics.contig_variable_overrides.json

# Build chromosome-specific override json of VCFs and VCF indexes
add_contig_vcfs_to_chromshard_overrides_json \
  $staging_dir/CollectVcfQcMetrics.contig_variable_overrides.json \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart2 \
  filtered_vcfs \
  filtered_vcf_idxs

# Write template input .json for short variant QC metric collection
cat << EOF > $staging_dir/CollectVcfQcMetrics.inputs.template.json
{
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": \$CONTIG_SCATTER_COUNT,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.001,
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:a5910a5",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.indel_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.indel.sites.bed.gz"],
  "CollectVcfQcMetrics.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 1000,
  "CollectVcfQcMetrics.output_prefix": "dfci-g2c.v1.gatkhc.initial_qc.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 7.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.snv.sites.bed.gz"],
  "CollectVcfQcMetrics.sv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.sv.sites.bed.gz"],
  "CollectVcfQcMetrics.trios_fam_file": "$MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam",
  "CollectVcfQcMetrics.vcfs": \$CONTIG_VCFS,
  "CollectVcfQcMetrics.vcf_idxs": \$CONTIG_VCF_IDXS
}
EOF

# Submit, monitor, stage, and cleanup short variant QC metadata workflow
code/scripts/manage_chromshards.py \
  --wdl code/wdl/pancan_germline_wgs/vcf-qc/CollectVcfQcMetrics.wdl \
  --input-json-template $staging_dir/CollectVcfQcMetrics.inputs.template.json \
  --contig-variable-overrides $staging_dir/CollectVcfQcMetrics.contig_variable_overrides.json \
  --dependencies-zip g2c.dependencies.zip \
  --staging-bucket $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/VcfQcMetrics/ \
  --name CollectInitialVcfQcMetrics \
  --contig-list contig_lists/dfci-g2c.v1.contigs.w$WN.list \
  --status-tsv cromshell/progress/dfci-g2c.v1.CollectVcfQcMetrics.initial_qc.progress.tsv \
  --workflow-id-log-prefix "dfci-g2c.v1" \
  --no-cleanup \
  --hard-reset \
  --outer-gate 30 \
  --max-attempts 1


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
  
  # Localize output tracker json and get URIs for QC metrics
  json_fname=CollectInitialVcfQcMetrics.chr$k.outputs.json
  gsutil cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/VcfQcMetrics/chr$k/$json_fname \
    $staging_dir/
  for key in size_distrib af_distrib size_vs_af_distrib \
             all_svs_bed common_snvs_bed common_indels_bed common_svs_bed; do
    jq .$key $staging_dir/$json_fname | fgrep -xv "null" | tr -d '"' >> $staging_dir/$key.uris.list
  done

  # Clear local copy of output tracker json
  rm $staging_dir/$json_fname
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
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:8009f0b",
  "PlotVcfQcMetrics.output_prefix": "dfci-g2c.v1.initial_qc",
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

# Once patches are complete, stage output
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.PlotInitialVcfQcMetrics.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/PlotQc/

# # Clear Cromwell execution & output buckets for patch jobs
# gsutil -m ls $( cat cromshell/job_ids/GnarlyJointGenotypingPart1.inputs.$contig.patch.job_ids.list \
#                 | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell/*/GnarlyJointGenotypingPart1/" \
#                   '{ print bucket_prefix$1"/**" }' ) \
# > uris_to_delete.list
# cleanup_garbage

