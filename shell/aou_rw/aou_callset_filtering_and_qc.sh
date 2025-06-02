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


###################################################
# Identify candidate twins / technical replicates #
###################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/fam_curation
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Build list of all autosomal GATK-HC VCFs
for suf in vcfs tbis; do
  flist=$staging_dir/gatkhc.$suf.uris.list
  if [ -e $flist ]; then rm $flist; fi
done
for k in $( seq 1 22 ); do
  contig="chr$k"

  json_fname="PosthocCleanupPart2.$contig.outputs.json"
  gsutil -m cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/PosthocCleanupPart2/$contig/$json_fname \
    $staging_dir/
  
  jq .filtered_vcfs $staging_dir/$json_fname \
  | sed 's/,/\n/g' | tr -d '[]"' | sed '/^$/d' | awk '{ print $1 }' \
  >> $staging_dir/gatkhc.vcfs.uris.list
  
  
  jq .filtered_vcf_idxs $staging_dir/$json_fname \
  | sed 's/,/\n/g' | tr -d '[]"' | sed '/^$/d' | awk '{ print $1 }' \
  >> $staging_dir/gatkhc.tbis.uris.list
done

# Write input .json for short variant QC metric collection
cat << EOF | python -m json.tool > cromshell/inputs/InferTwins.inputs.json
{
  "InferTwins.ConcatVcfs.bcftools_concat_options": "--allow-overlaps",
  "InferTwins.ConcatVcfs.cpu_cores": 4,
  "InferTwins.ConcatVcfs.disk_gb": 1000,
  "InferTwins.ConcatVcfs.mem_gb": 16,
  "InferTwins.PrepVcf.mem_gb": 8.0,
  "InferTwins.PrepVcf.n_cpu": 4,
  "InferTwins.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "InferTwins.min_kin": 0.45,
  "InferTwins.only_snvs": true,
  "InferTwins.out_prefix": "dfci-g2c.v1",
  "InferTwins.plink2_docker": "vanallenlab/g2c_pipeline:plink2",
  "InferTwins.vcfs": $( collapse_txt $staging_dir/gatkhc.vcfs.uris.list ),
  "InferTwins.vcf_idxs": $( collapse_txt $staging_dir/gatkhc.tbis.uris.list )
}
EOF

# Submit twin inference workflow
cromshell --no_turtle -t 120 -mc submit \
  --options-json code/refs/json/aou.cromwell_options.default.json \
  --dependencies-zip g2c.dependencies.zip \
  code/wdl/pancan_germline_wgs/InferTwins.wdl \
  cromshell/inputs/InferTwins.inputs.json \
| jq .id | tr -d '"' \
>> cromshell/job_ids/dfci-g2c.v1.InferTwins.job_ids.list

# Monitor QC visualization workflow
monitor_workflow $( tail -n1 cromshell/job_ids/dfci-g2c.v1.InferTwins.job_ids.list )

# Once workflow is complete, stage output
cromshell -t 120 list-outputs \
  $( tail -n1 cromshell/job_ids/dfci-g2c.v1.InferTwins.job_ids.list ) \
| awk '{ print $2 }' \
| gsutil -m cp -I \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/

# Clear Cromwell execution & output buckets
gsutil -m ls $( cat cromshell/job_ids/dfci-g2c.v1.InferTwins.job_ids.list \
                | awk -v bucket_prefix="$WORKSPACE_BUCKET/cromwell-*/InferTwins/" \
                  '{ print bucket_prefix$1"/**" }' ) \
> uris_to_delete.list
cleanup_garbage

# Locally filter putative twins to remove pairs that are more likely to be 
# duplicated WGS files (based on WGS QC metric similarity) than true replicates
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.kin0.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/
code/scripts/clean_candidate_twins.R \
  $staging_dir/dfci-g2c.v1.kin0.gz \
  $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/dfci-g2c.v1.cleaned.kin0
gzip -f $staging_dir/dfci-g2c.v1.cleaned.kin0
gsutil -m cp \
  $staging_dir/dfci-g2c.v1.cleaned.kin0.gz \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/


#########################################################################
# Curate existing long-read and short-read WGS callsets for AoU samples #
#########################################################################

# TODO: implement this



###############################################################
# Build QC sample inclusion priority & sampling probabilities #
###############################################################

# Note: this only needs to be run once for the entire cohort across all workspaces

# Reaffirm staging directory
staging_dir=staging/sample_priority
if ! [ -e $staging_dir ]; then mkdir $staging_dir; fi

# Localize most up-to-date sample QC manifest
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-hc/qc-filtering/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  $staging_dir/

# Generate sampling probabilities
code/scripts/assign_sample_qc_weights.R \
  --qc-tsv $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
  --pass-column global_qc_pass \
  --pass-column batch_qc_pass \
  --pass-column clusterbatch_qc_pass \
  --pass-column filtersites_qc_pass \
  --pass-column gatksv_posthoc_qc_pass \
  --pass-column gatkhc_posthoc_qc_pass \
  --out-tsv $staging_dir/dfci-g2c.v1.qc.sample_weights.tsv

# Extract list of sample IDs still present in callset at this stage
zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| awk -v FS="\t" -v OFS="\t" '{ if ($NF=="True") print $1 }' \
> $staging_dir/g2c_ids.present_after_calling.samples.list
zcat $staging_dir/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz \
| awk -v FS="\t" -v OFS="\t" '{ if ($NF=="True") print $1, $2, $3 }' \
> $staging_dir/g2c_ids.present_after_calling.sample_cohort_map.tsv

# Find list of samples included in complete trios
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/data/sample_info/relatedness/dfci-g2c.reported_families.fam \
  $staging_dir/
code/scripts/subset_fam.R \
  --in-fam $staging_dir/dfci-g2c.reported_families.fam \
  --all-samples $staging_dir/g2c_ids.present_after_calling.samples.list \
  --out-fam $staging_dir/dfci-g2c.reported_families.filtered.fam
awk -v FS="\t" -v OFS="\t" '{ if ($2!=0 && $3!=0 && $4!=0) print }' \
  $staging_dir/dfci-g2c.reported_families.filtered.fam \
> $staging_dir/dfci-g2c.reported_families.filtered.complete.fam
awk -v OFS="\n" '{ print $2, $3, $4 }' \
  $staging_dir/dfci-g2c.reported_families.filtered.complete.fam \
| sort -V | uniq \
> $staging_dir/dfci-g2c.reported_families.filtered.complete.samples.list

# Find list of samples included in complete twin pairs
gsutil -m cp \
  $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/InferTwins/dfci-g2c.v1.cleaned.kin0.gz \
  $staging_dir/
zcat $staging_dir/dfci-g2c.v1.cleaned.kin0.gz \
| fgrep -v "#" | awk -v OFS="\t" '{ print "twins_"NR, $2, $4, 0, 0 }' \
> $staging_dir/all_twins.fam
code/scripts/subset_fam.R \
  --in-fam $staging_dir/all_twins.fam \
  --all-samples $staging_dir/g2c_ids.present_after_calling.samples.list \
  --out-fam $staging_dir/all_twins.filtered.fam
awk -v FS="\t" -v OFS="\n" '{ if ($2!=0 && $3!=0) print $2, $3 }' \
  $staging_dir/all_twins.filtered.fam \
| sort -V | uniq \
> $staging_dir/dfci-g2c.inferred_twins.complete.samples.list

# Get list of samples with complete HGSV long-read WGS calls
# Supplement with twins/replicates, since these are equally useful for our purposes
gsutil -m cat \
  gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/snv_indel/1KGP.lrWGS.snv_indel.cleaned.chrY.vcf.gz \
| gunzip -c | head -n10000 | fgrep "#" | bcftools query -l \
> $staging_dir/1KGP.lrWGS.snv_indel.external_ids.tsv
gsutil -m cat \
  gs://dfci-g2c-refs/hgsv/dense_vcfs/lrwgs/sv/1KGP.lrWGS.sv.cleaned.chrY.vcf.gz \
| gunzip -c | head -n10000 | fgrep "#" | bcftools query -l \
> $staging_dir/1KGP.lrWGS.sv.external_ids.tsv
fgrep -xf \
  $staging_dir/1KGP.lrWGS.snv_indel.external_ids.tsv \
  $staging_dir/1KGP.lrWGS.sv.external_ids.tsv \
| sort | uniq \
> $staging_dir/1KGP.lrWGS.complete.external_ids.tsv
awk -v OFS="\t" '{ if ($3=="hgsvc") print $2, $1 }' \
  $staging_dir/g2c_ids.present_after_calling.sample_cohort_map.tsv \
| sort -k1,1 \
| join -j 1 -t $'\t' - $staging_dir/1KGP.lrWGS.complete.external_ids.tsv \
| cut -f2 | sort -V \
> $staging_dir/1KGP.lrWGS.complete.g2c_ids.no_twins.samples.list
zcat $staging_dir/dfci-g2c.v1.cleaned.kin0.gz | cut -f2,4 \
| fgrep -wf $staging_dir/1KGP.lrWGS.complete.g2c_ids.no_twins.samples.list \
| sed 's/\t/\n/g' \
| cat - $staging_dir/1KGP.lrWGS.complete.g2c_ids.no_twins.samples.list \
| sort -V | uniq \
> $staging_dir/1KGP.lrWGS.complete.g2c_ids.samples.list

# Get list of samples with complete AoU long-read WGS calls
# Supplement with twins/replicates, since these are equally useful for our purposes
# TODO: implement this

# Combine lists of samples with lrWGS available from at least one source
# TODO: implement this
cp \
  $staging_dir/1KGP.lrWGS.complete.g2c_ids.samples.list \
  $staging_dir/dfci-g2c.lrWGS.complete.samples.list

# Compute sample priority and append sampling probabilities
cat \
  $staging_dir/dfci-g2c.reported_families.filtered.complete.samples.list \
  $staging_dir/dfci-g2c.inferred_twins.complete.samples.list \
  $staging_dir/dfci-g2c.lrWGS.complete.samples.list \
  $staging_dir/g2c_ids.present_after_calling.samples.list \
| sort -V | uniq -c | awk -v OFS="\t" '{ print $2, $1 }' | sort -k1,1 \
| join -j 1 -t $'\t' \
  - <( sort -k1,1 $staging_dir/dfci-g2c.v1.qc.sample_weights.tsv ) \
| sort -nrk2,2 -k3,3nr -k1,1V \
| cat <( echo -e "#G2C_ID\tpriority\tweight" ) - \
> $staging_dir/dfci-g2c.v1.sample_qc_priority.tsv


#####################################
# Build inter-cohort sample ID maps #
#####################################

# TODO: implement this

# # Build map of sample IDs for 1KGP srWGS short variants
# gsutil -m cat \
#   gs://dfci-g2c-refs/hgsv/dense_vcfs/srwgs/snv_indel/1KGP.srWGS.snv_indel.cleaned.chrY.vcf.gz \
# | gunzip -c | head -n10000 | fgrep "#" | bcftools query -l \
# > $staging_dir/1KGP.srWGS.snv_indel.external_ids.tsv
# # TODO: write script to make alias table, accounting for putative twins


##############################
# Collect initial QC metrics #
##############################

# Note: this workflow below is scattered across all five workspaces for 
# max parallelization. It must be submitted as below in each workspace.

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

# Write template input .json for QC metric collection
cat << EOF > $staging_dir/CollectVcfQcMetrics.inputs.template.json
{
  "CollectVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "CollectVcfQcMetrics.benchmarking_shards": \$CONTIG_SCATTER_COUNT,
  "CollectVcfQcMetrics.benchmark_interval_beds": ["gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.easy.\$CONTIG.bed.gz",
                                                  "gs://dfci-g2c-refs/giab/\$CONTIG/giab.hg38.broad_callable.hard.\$CONTIG.bed.gz"],
  "CollectVcfQcMetrics.benchmark_interval_bed_names": ["giab_easy", "giab_hard"],
  "CollectVcfQcMetrics.common_af_cutoff": 0.001,
  "CollectVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:d8cdae4",
  "CollectVcfQcMetrics.genome_file": "gs://dfci-g2c-refs/hg38/hg38.genome",
  "CollectVcfQcMetrics.linux_docker": "marketplace.gcr.io/google/ubuntu1804",
  "CollectVcfQcMetrics.n_for_sample_level_analyses": 1000,
  "CollectVcfQcMetrics.output_prefix": "dfci-g2c.v1.initial_qc.\$CONTIG",
  "CollectVcfQcMetrics.PreprocessVcf.mem_gb": 7.5,
  "CollectVcfQcMetrics.PreprocessVcf.n_cpu": 4,
  "CollectVcfQcMetrics.shard_vcf": false,
  "CollectVcfQcMetrics.site_benchmark_dataset_names": ["gnomad_v4"],
  "CollectVcfQcMetrics.snv_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.snv.sites.bed.gz"],
  "CollectVcfQcMetrics.indel_site_benchmark_beds": ["gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/\$CONTIG/gnomad.v4.1.\$CONTIG.indel.sites.bed.gz"],
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
  fname=$staging_dir/$key.uris.list
  if [ -e $fname ]; then rm $fname; fi
done
for key in site_benchmark_ppv_by_freqs site_benchmark_sensitivity_by_freqs \
           site_benchmark_common_snv_ppv_beds \
           site_benchmark_common_indel_ppv_beds \
           site_benchmark_common_sv_ppv_beds; do
  for subset in giab_easy giab_hard; do
    fname=$staging_dir/$key.$subset.uris.list
    if [ -e $fname ]; then rm $fname; fi
  done
done
for suffix in af_distribution size_distribution; do
  fname=$staging_dir/gnomAD_$suffix.list
  if [ -e $fname ]; then rm $fname; fi
done
for k in $( seq 1 22 ) X Y; do
  
  # Localize output tracker json and get URIs for QC metrics
  json_fname=CollectInitialVcfQcMetrics.chr$k.outputs.json
  gsutil cp \
    $MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/qc-filtering/initial-qc/VcfQcMetrics/chr$k/$json_fname \
    $staging_dir/
  for key in size_distrib af_distrib size_vs_af_distrib \
             all_svs_bed common_snvs_bed common_indels_bed common_svs_bed; do
    jq .\"CollectVcfQcMetrics.$key\" $staging_dir/$json_fname \
    | fgrep -xv "null" | tr -d '"' \
    >> $staging_dir/$key.uris.list
  done
  for key in site_benchmark_ppv_by_freqs site_benchmark_sensitivity_by_freqs \
             site_benchmark_common_snv_ppv_beds \
             site_benchmark_common_indel_ppv_beds \
             site_benchmark_common_sv_ppv_beds; do
    for subset in giab_easy giab_hard; do
      jq .\"CollectVcfQcMetrics.$key\" $staging_dir/$json_fname \
      | fgrep -xv "null" | tr -d '"[]' | sed 's/,$/\n/g' \
      | sed '/^$/d' | awk '{ print $1 }' | fgrep $subset \
      >> $staging_dir/$key.$subset.uris.list
    done
  done

  # Clear local copy of output tracker json
  rm $staging_dir/$json_fname

  # Add precomputed gnomAD v4.1 reference distributions to file lists
  for suffix in af_distribution size_distribution; do
    echo "gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/chr$k/gnomad.v4.1.chr$k.$suffix.merged.tsv.gz" \
    >> $staging_dir/gnomAD_$suffix.uris.list
  done
done

# Write input .json for short variant QC metric collection
cat << EOF | python -m json.tool > cromshell/inputs/PlotInitialVcfQcMetrics.inputs.json
{
  "PlotVcfQcMetrics.af_distribution_tsvs": $( collapse_txt $staging_dir/af_distrib.uris.list ),
  "PlotVcfQcMetrics.all_sv_beds": $( collapse_txt $staging_dir/all_svs_bed.uris.list ),
  "PlotVcfQcMetrics.bcftools_docker": "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base-mini:2024-10-25-v0.29-beta-5ea22a52",
  "PlotVcfQcMetrics.common_af_cutoff": 0.001,
  "PlotVcfQcMetrics.common_snv_beds": $( collapse_txt $staging_dir/common_snvs_bed.uris.list ),
  "PlotVcfQcMetrics.common_indel_beds": $( collapse_txt $staging_dir/common_indels_bed.uris.list ),
  "PlotVcfQcMetrics.common_sv_beds": $( collapse_txt $staging_dir/common_svs_bed.uris.list ),
  "PlotVcfQcMetrics.g2c_analysis_docker": "vanallenlab/g2c_analysis:b9557dc",
  "PlotVcfQcMetrics.output_prefix": "dfci-g2c.v1.initial_qc",
  "PlotVcfQcMetrics.ref_af_distribution_tsvs": $( collapse_txt $staging_dir/gnomAD_af_distribution.uris.list ),
  "PlotVcfQcMetrics.ref_size_distribution_tsvs": $( collapse_txt $staging_dir/gnomAD_size_distribution.uris.list ),
  "PlotVcfQcMetrics.ref_cohort_prefix": "gnomAD_v4.1",
  "PlotVcfQcMetrics.ref_cohort_plot_title": "gnomAD v4.1",
  "PlotVcfQcMetrics.site_benchmark_common_snv_ppv_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_snv_ppv_beds.giab_easy.uris.list ),
                                                            $( collapse_txt $staging_dir/site_benchmark_common_snv_ppv_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_common_indel_ppv_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_indel_ppv_beds.giab_easy.uris.list ),
                                                              $( collapse_txt $staging_dir/site_benchmark_common_indel_ppv_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_common_sv_ppv_beds": [[ $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_easy.uris.list ),
                                                           $( collapse_txt $staging_dir/site_benchmark_common_sv_ppv_beds.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_ppv_by_freqs": [[ $( collapse_txt $staging_dir/site_benchmark_ppv_by_freqs.giab_easy.uris.list ),
                                                     $( collapse_txt $staging_dir/site_benchmark_ppv_by_freqs.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_sensitivity_by_freqs": [[ $( collapse_txt $staging_dir/site_benchmark_sensitivity_by_freqs.giab_easy.uris.list ),
                                                             $( collapse_txt $staging_dir/site_benchmark_sensitivity_by_freqs.giab_hard.uris.list ) ]],
  "PlotVcfQcMetrics.site_benchmark_dataset_prefixes": ["gnomad_v4.1"],
  "PlotVcfQcMetrics.site_benchmark_dataset_titles": ["gnomAD v4.1"],
  "PlotVcfQcMetrics.site_benchmark_interval_names": ["Easy", "Hard"],
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

# Once workflow is complete, stage output
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

