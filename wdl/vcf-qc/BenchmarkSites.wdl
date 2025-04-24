# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform two-sided site-level benchmarking of a "source" callset versus a reference "target" callset

# Expects that both callsets have already been processed by CollectVcfQcMetrics.wdl


version 1.0


import "BenchmarkSitesSingle.wdl" as BenchSingle
import "Utilities.wdl" as Utils


workflow BenchmarkSites {
  input {
    File? source_snv_bed
    File? source_indel_bed
    File? source_sv_bed
    String source_prefix

    File? target_snv_bed
    File? target_indel_bed
    File? target_sv_bed
    String target_prefix

    Array[File] eval_interval_beds
    Array[String]? eval_interval_bed_names
    File genome_file

    Int total_shards = 2500
    Float common_af_cutoff

    String bcftools_docker
    String g2c_analysis_docker
  }

  # Index SNV BEDs
  if ( defined(source_snv_bed) ){
    call Utils.IndexBed as IndexSourceSnvs {
      input:
        bed = select_first([source_snv_bed]),
        docker = bcftools_docker
    }
  }
  if ( defined(target_snv_bed) ){
    call Utils.IndexBed as IndexTargetSnvs {
      input:
        bed = select_first([target_snv_bed]),
        docker = bcftools_docker
    }
  }

  # Index indel BEDs
  if ( defined(source_indel_bed) ){
    call Utils.IndexBed as IndexSourceIndels {
      input:
        bed = select_first([source_indel_bed]),
        docker = bcftools_docker
    }
  }
  if ( defined(target_indel_bed) ){
    call Utils.IndexBed as IndexTargetIndels {
      input:
        bed = select_first([target_indel_bed]),
        docker = bcftools_docker
    }
  }

  # Index SV BEDs
  if ( defined(source_sv_bed) ){
    call Utils.IndexBed as IndexSourceSvs {
      input:
        bed = select_first([source_sv_bed]),
        docker = bcftools_docker
    }
  }
  if ( defined(target_sv_bed) ){
    call Utils.IndexBed as IndexTargetSvs {
      input:
        bed = select_first([target_sv_bed]),
        docker = bcftools_docker
    }
  }

  Int n_eval_beds = length(eval_interval_beds)
  Int shards_per_eval_bed = floor(total_shards / n_eval_beds)

  scatter ( eval_idx in range(n_eval_beds) ){

    File eval_bed = eval_interval_beds[eval_idx]
    if ( defined(eval_interval_bed_names) ) {
      String eval_name_custom = select_first([eval_interval_bed_names])[eval_idx]
    }
    String eval_name = select_first([eval_name_custom, basename(eval_bed, ".bed.gz")])

    # Run evaluation for a single eval BED as a subworkflow
    call BenchSingle.BenchmarkSitesSingle as BenchmarkTask {
      input:
        eval_interval_bed = eval_bed,
        genome_file = genome_file,
        source_snv_bed = source_snv_bed,
        source_snv_bed_idx = IndexSourceSnvs.bed_idx,
        target_snv_bed = target_snv_bed,
        target_snv_bed_idx = IndexTargetSnvs.bed_idx,
        source_indel_bed = source_indel_bed,
        source_indel_bed_idx = IndexSourceIndels.bed_idx,
        target_indel_bed = target_indel_bed,
        target_indel_bed_idx = IndexTargetIndels.bed_idx,
        source_sv_bed = source_sv_bed,
        source_sv_bed_idx = IndexSourceSvs.bed_idx,
        target_sv_bed = target_sv_bed,
        target_sv_bed_idx = IndexTargetSvs.bed_idx,
        eval_prefix = eval_name,
        source_prefix = source_prefix,
        target_prefix = target_prefix,
        common_af_cutoff = common_af_cutoff,
        total_shards = shards_per_eval_bed,
        bcftools_docker = bcftools_docker,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  output {
    Array[File?] common_snv_ppv_beds = BenchmarkTask.common_snv_ppv_bed
    Array[File?] common_snv_sens_beds = BenchmarkTask.common_indel_ppv_bed
    Array[File?] common_indel_ppv_beds = BenchmarkTask.common_indel_ppv_bed
    Array[File?] common_indel_sens_beds = BenchmarkTask.common_indel_sens_bed
    Array[File?] common_sv_ppv_beds = BenchmarkTask.common_sv_ppv_bed
    Array[File?] common_sv_sens_beds = BenchmarkTask.common_sv_sens_bed
    Array[File] ppv_by_sizes = BenchmarkTask.ppv_by_size
    Array[File] sensitivity_by_sizes = BenchmarkTask.sensitivity_by_size
    Array[File] ppv_by_freqs = BenchmarkTask.ppv_by_freq
    Array[File] sensitivity_by_freqs = BenchmarkTask.sensitivity_by_freq
  }
}

