import os
import subprocess
import json

sample_name_map="gs://fc-secure-d531c052-7b41-4dea-9e1d-22e648f6e228/ufc_jg/jg_part1/dfci-ufc.sample_map"

# Construct the base path for the Google Cloud Storage bucket
bucket_path = os.getenv('WORKSPACE_BUCKET')

# Function to list paths using gsutil and construct arrays
def process_directory(dir_path,suffix):
    command = f"gsutil ls {bucket_path}/GnarlyJointGenotypingPart1-Output/{callset_name}/{dir_path}/*{suffix}"
    paths = subprocess.check_output(command, shell=True).decode('utf-8').strip().split('\n')
    return paths

# List all directories under ${bucket_path}/GnarlyJointGenotypingPart1-Output/
command = f"gsutil ls -d {bucket_path}/GnarlyJointGenotypingPart1-Output/*/"
callset_names = subprocess.check_output(command, shell=True).decode('utf-8').strip().split('\n')
callset_names = [name.split('/')[-2] for name in callset_names]

# Instantiate Arrays
output_sites_only_vcf = []
output_variant_filtered_vcf = []
output_sites_only_vcf_idx = []
output_variant_filtered_vcf_idx = []

# Loop through each callset_name
for callset_name in callset_names:
    # Process each directory and populate respective arrays
    sites_only_vcf = process_directory("sites_only_vcf",".vcf.gz")
    variant_filtered_vcf = process_directory("variant_filtered_vcf",".vcf.gz")

    # Process each directory and populate respective arrays
    sites_only_vcf_idx = process_directory("sites_only_vcf",".tbi")
    variant_filtered_vcf_idx = process_directory("variant_filtered_vcf",".tbi")

    # Append paths to respective output lists
    output_sites_only_vcf.append(sites_only_vcf)
    output_variant_filtered_vcf.append(variant_filtered_vcf)
    output_sites_only_vcf_idx.append(sites_only_vcf_idx)
    output_variant_filtered_vcf_idx.append(variant_filtered_vcf_idx)
# JSON file path
output_json_file = "GnarlyJointGenotypingPart2.json"

# Construct JSON content
json_content = {
    "GnarlyJointGenotypingPart2.ApplyRecalibration.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.CollectMetricsOnFullVcf.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.CollectMetricsSharded.gatk-publ_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.FinalGatherVcf.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.GatherVariantCallingMetrics.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.IndelsVariantRecalibrator.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.SNPGatherTranches.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.SNPsVariantRecalibratorCreateModel.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.SNPsVariantRecalibratorScattered.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.SitesOnlyGatherVcf.gatk_docker": "broadinstitute/gatk:4.4.0.0",
    "GnarlyJointGenotypingPart2.axiomPoly_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz",
    "GnarlyJointGenotypingPart2.axiomPoly_resource_vcf_index": "gs://gcp-public-data--broad-references/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi",
    "GnarlyJointGenotypingPart2.callset_name": "dfci-ufc",
    "GnarlyJointGenotypingPart2.dbsnp_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf",
    "GnarlyJointGenotypingPart2.dbsnp_resource_vcf_index": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
    "GnarlyJointGenotypingPart2.eval_interval_list": "gs://gcp-public-data--broad-references/hg38/v0/wgs_evaluation_regions.hg38.interval_list",
    "GnarlyJointGenotypingPart2.hapmap_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz",
    "GnarlyJointGenotypingPart2.hapmap_resource_vcf_index": "gs://gcp-public-data--broad-references/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi",
    "GnarlyJointGenotypingPart2.mills_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz",
    "GnarlyJointGenotypingPart2.mills_resource_vcf_index": "gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi",
    "GnarlyJointGenotypingPart2.omni_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz",
    "GnarlyJointGenotypingPart2.omni_resource_vcf_index": "gs://gcp-public-data--broad-references/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi",
    "GnarlyJointGenotypingPart2.one_thousand_genomes_resource_vcf": "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz",
    "GnarlyJointGenotypingPart2.one_thousand_genomes_resource_vcf_index": "gs://gcp-public-data--broad-references/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi",
    "GnarlyJointGenotypingPart2.ref_dict": "gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict",
    "GnarlyJointGenotypingPart2.sites_only_vcfs_index_nested": output_sites_only_vcf_idx,
    "GnarlyJointGenotypingPart2.sites_only_vcfs_nested": output_sites_only_vcf,
    "GnarlyJointGenotypingPart2.variant_filtered_vcfs_index_nested": output_variant_filtered_vcf_idx,
    "GnarlyJointGenotypingPart2.variant_filtered_vcfs_nested": output_variant_filtered_vcf,
    "GnarlyJointGenotypingPart2.sample_name_map": sample_name_map
}

# Write JSON content to file without prettifying
with open(output_json_file, 'w') as f:
    json.dump(json_content, f, separators=(',', ':'))

    # Add a newline character at the end of the file
    f.write('\n')

print(f"Output written to {output_json_file}")