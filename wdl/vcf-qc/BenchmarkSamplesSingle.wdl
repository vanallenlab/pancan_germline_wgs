# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform two-sided sample-level benchmarking of a "source" callset versus a reference "target" callset
# Note: this is a subroutine called within BenchmarkSamples.wdl, and is designed to be
# invoked by this parent workflow

# Expects that both callsets are provided according to G2C QC expectations


version 1.0


import "BenchmarkSites.wdl" as SiteBench
import "QcTasks.wdl" as QcTasks


workflow BenchmarkSamplesSingle {
  input {
    File? source_snv_bed
    File? source_indel_bed
    File? source_sv_bed
    File source_gt_tarball
    File source_sample_id_list
    String source_prefix

    File target_vcf
    File target_vcf_idx
    String target_prefix
    
    File id_map_tsv

    Array[File] eval_interval_beds
    Array[String]? eval_interval_bed_names
    File genome_file

    Int total_shards = 2500
    Float common_af_cutoff = 0.01

    String bcftools_docker
    String g2c_analysis_docker
  }

  # Subset target VCF to overlapping samples
  call SubsetTargetVcf {
    input:
      vcf = target_vcf,
      vcf_idx = target_vcf_idx,
      id_map_tsv = id_map_tsv,
      source_sample_id_list = source_sample_id_list,
      bcftools_docker = bcftools_docker
  }

  # Collect site metrics for target VCF
  call QcTasks.CollectSiteMetrics as TargetVcfToBed {
    input:
      vcf = SubsetTargetVcf.filtered_sites_vcf,
      vcf_idx = SubsetTargetVcf.filtered_sites_vcf_idx,
      n_samples = SubsetTargetVcf.n_samples_in_original_vcf,
      common_af_cutoff = common_af_cutoff,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Run site-level benchmarking to determine mapping of 
  # matching variant IDs between source & target
  call SiteBench.BenchmarkSites {
    input:
      source_snv_bed = source_snv_bed,
      source_indel_bed = source_indel_bed,
      source_sv_bed = source_sv_bed,
      source_prefix = source_prefix,
      target_snv_bed = TargetVcfToBed.snv_sites,
      target_indel_bed = TargetVcfToBed.indel_sites,
      target_sv_bed = TargetVcfToBed.sv_sites,
      target_prefix = target_prefix,
      eval_interval_beds = eval_interval_beds,
      eval_interval_bed_names = eval_interval_bed_names,
      genome_file = genome_file,
      total_shards = total_shards,
      common_af_cutoff = common_af_cutoff,
      make_full_id_maps = true,
      bcftools_docker = bcftools_docker,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collect sample genotypes for subsetted target VCF
  call QcTasks.CollectSampleGenotypeMetrics as GatherTargetGts {
    input:
      vcf = SubsetTargetVcf.filtered_dense_vcf,
      vcf_idx = SubsetTargetVcf.filtered_dense_vcf_idx,
      site_metrics = TargetVcfToBed.all_sites,
      g2c_analysis_docker = g2c_analysis_docker
  }
  call QcTasks.ConcatGenotypeTsvs as MakeTargetGtTarball {
    input:
      tsvs = [GatherTargetGts.genotypes_tsv],
      output_prefix = target_prefix,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # TODO:
  # 5. Compare .tsvs from 4 to equivalent .tsvs generated from target/input callset while linking through variant ID map from #3
  # 6. Collapse outputs of 5 across all samples and compute summary metrics (medians?)
  
  output {}
}


# Subset a target VCF based on a dynamically-determined list of overlapping samples
task SubsetTargetVcf {
  input {
    File vcf
    File vcf_idx
    File id_map_tsv
    File source_sample_id_list
    String bcftools_docker
  }

  Int disk_gb = ceil(2.5 * size(vcf, "GB")) + 10

  String filtered_map_fname = basename(id_map_tsv, ".tsv") + ".filtered.tsv"
  String filtered_dense_vcf_fname = basename(vcf, ".vcf.gz") + "filtered.vcf.gz"
  String filtered_sites_vcf_fname = basename(vcf, ".vcf.gz") + "filtered.sites.vcf.gz"

  command <<<
    set -eu -o pipefail

    # Get list of samples present in input VCF
    bcftools query -l ~{vcf} | sort > target.samples.list

    # Ensure all files are sorted for bash join
    sort ~{source_sample_id_list} > source.samples.list
    sort -k1,1 ~{id_map_tsv} > id_map.tsv

    # Use bash join to subset ID map
    join -j 1 -t $'\t' source.samples.list id_map.tsv \
    | sort -k2,2 \
    | join -1 2 -2 1 -t $'\t' - target.samples.list \
    > ~{filtered_map_fname}

    # Subset input VCF to overlapping samples
    cut -f2 ~{filtered_map_fname} | sort -V | uniq > target.keep.list
    bcftools view \
      --no-update \
      --samples-file target.keep.list \
      --force-samples \
      ~{vcf} \
    | bcftools view \
      --no-update \
      --include 'GT="alt" | FILTER="MULTIALLELIC"' \
      -Oz -o ~{filtered_dense_vcf_fname}
    tabix -p vcf -f ~{filtered_dense_vcf_fname}

    # Make a sites VCF for benchmarking comparison
    bcftools view --no-update -G \
      -Oz -o ~{filtered_sites_vcf_fname} \
      ~{filtered_dense_vcf_fname}
    tabix -p vcf -f ~{filtered_sites_vcf_fname}
  >>>

  output {
    File filtered_dense_vcf = filtered_dense_vcf_fname
    File filtered_dense_vcf_idx = "~{filtered_dense_vcf_fname}.tbi"

    File filtered_sites_vcf = filtered_sites_vcf_fname
    File filtered_sites_vcf_idx = "~{filtered_sites_vcf_fname}.tbi"

    File filtered_id_map = filtered_map_fname

    Int n_samples_in_original_vcf = length(read_lines("target.samples.list"))
  }

  runtime {
    docker: bcftools_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}

