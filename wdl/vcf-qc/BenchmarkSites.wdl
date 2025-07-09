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

    Boolean make_full_id_maps = false

    String bcftools_docker
    String g2c_analysis_docker
  }


  # Index SNV BEDs
  if ( defined(source_snv_bed) ){
    call Utils.MakeTabixIndex as IndexSourceSnvs {
      input:
        input_file = select_first([source_snv_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }
  if ( defined(target_snv_bed) ){
    call Utils.MakeTabixIndex as IndexTargetSnvs {
      input:
        input_file = select_first([target_snv_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }


  # Index indel BEDs
  if ( defined(source_indel_bed) ){
    call Utils.MakeTabixIndex as IndexSourceIndels {
      input:
        input_file = select_first([source_indel_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }
  if ( defined(target_indel_bed) ){
    call Utils.MakeTabixIndex as IndexTargetIndels {
      input:
        input_file = select_first([target_indel_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }


  # Index SV BEDs
  if ( defined(source_sv_bed) ){
    call Utils.MakeTabixIndex as IndexSourceSvs {
      input:
        input_file = select_first([source_sv_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }
  if ( defined(target_sv_bed) ){
    call Utils.MakeTabixIndex as IndexTargetSvs {
      input:
        input_file = select_first([target_sv_bed]),
        file_type = "bed",
        docker = bcftools_docker
    }
  }


  # Main benchmarking block
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
        source_snv_bed_idx = IndexSourceSnvs.tbi,
        target_snv_bed = target_snv_bed,
        target_snv_bed_idx = IndexTargetSnvs.tbi,
        source_indel_bed = source_indel_bed,
        source_indel_bed_idx = IndexSourceIndels.tbi,
        target_indel_bed = target_indel_bed,
        target_indel_bed_idx = IndexTargetIndels.tbi,
        source_sv_bed = source_sv_bed,
        source_sv_bed_idx = IndexSourceSvs.tbi,
        target_sv_bed = target_sv_bed,
        target_sv_bed_idx = IndexTargetSvs.tbi,
        eval_prefix = eval_name,
        source_prefix = source_prefix,
        target_prefix = target_prefix,
        common_af_cutoff = common_af_cutoff,
        total_shards = shards_per_eval_bed,
        make_full_id_maps = make_full_id_maps,
        bcftools_docker = bcftools_docker,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }


  # Perform false positive run detection
  Array[File] all_common_ppv_beds = select_all(flatten([select_all(BenchmarkTask.common_snv_ppv_bed),
                                                        select_all(BenchmarkTask.common_indel_ppv_bed),
                                                        select_all(BenchmarkTask.common_sv_ppv_bed)]))

  call Utils.ConcatTextFiles as CollapseCommonPpvBeds {
    input:
      shards = all_common_ppv_beds,
      concat_command = "zcat",
      sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
      compression_command = "bgzip -c",
      input_has_header = true,
      output_filename = "all_common.ppv.bed.gz",
      docker = bcftools_docker
  }
  
  call DetectBadRuns as DetectFPRuns {
    input:
      bed = CollapseCommonPpvBeds.merged_file,
      output_prefix = "~{source_prefix}.~{target_prefix}.false_positive_runs",
      g2c_analysis_docker = g2c_analysis_docker
  }


  # Perform false negative run detection
  Array[File] all_common_sens_beds = select_all(flatten([select_all(BenchmarkTask.common_snv_sens_bed),
                                                        select_all(BenchmarkTask.common_indel_sens_bed),
                                                        select_all(BenchmarkTask.common_sv_sens_bed)]))

  call Utils.ConcatTextFiles as CollapseCommonSensBeds {
    input:
      shards = all_common_sens_beds,
      concat_command = "zcat",
      sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
      compression_command = "bgzip -c",
      input_has_header = true,
      output_filename = "all_common.sens.bed.gz",
      docker = bcftools_docker
  }

  call DetectBadRuns as DetectFNRuns {
    input:
      bed = CollapseCommonSensBeds.merged_file,
      output_prefix = "~{source_prefix}.~{target_prefix}.false_negative_runs",
      g2c_analysis_docker = g2c_analysis_docker
  }


  output {
    Array[File] common_snv_ppv_beds = select_all(BenchmarkTask.common_snv_ppv_bed)
    Array[File] common_snv_sens_beds = select_all(BenchmarkTask.common_indel_ppv_bed)
    Array[File] common_indel_ppv_beds = select_all(BenchmarkTask.common_indel_ppv_bed)
    Array[File] common_indel_sens_beds = select_all(BenchmarkTask.common_indel_sens_bed)
    Array[File] common_sv_ppv_beds = select_all(BenchmarkTask.common_sv_ppv_bed)
    Array[File] common_sv_sens_beds = select_all(BenchmarkTask.common_sv_sens_bed)
    Array[File] ppv_by_sizes = BenchmarkTask.ppv_by_size
    Array[File] sensitivity_by_sizes = BenchmarkTask.sensitivity_by_size
    Array[File] ppv_by_freqs = BenchmarkTask.ppv_by_freq
    Array[File] sensitivity_by_freqs = BenchmarkTask.sensitivity_by_freq
    Array[File?] ppv_variant_id_maps = BenchmarkTask.ppv_variant_id_map
    Array[File?] sensitivity_variant_id_maps = BenchmarkTask.sensitivity_variant_id_map
    File fp_runs = DetectFPRuns.runs_bed
    File fn_runs = DetectFNRuns.runs_bed
  }
}


task DetectBadRuns {
  input {
    File bed
    String output_prefix

    Int min_variants = 10
    Int min_run_length = 50000

    String g2c_analysis_docker
  }

  Int disk_gb = ceil(2 * size(bed, "GB")) + 10
  String out_fname = "~{output_prefix}.bed.gz"

  command <<<
    set -eu -o pipefail

    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/find_problematic_intervals.py \
      --min-vars ~{min_variants} \
      --min-bp ~{min_run_length} \
      ~{bed} \
      ~{out_fname}
  >>>

  output {
    File runs_bed = out_fname
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.7 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 2
    max_retries: 1
  }
}

