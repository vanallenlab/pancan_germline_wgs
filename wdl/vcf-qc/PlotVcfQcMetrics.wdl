# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot quality metrics for assessing a GATK joint-genotyped VCF (SNVs/indels or SVs)
# See CollectVcfQcMetrics.wdl for the generation of inputs for this workflow


version 1.0


import "QcTasks.wdl" as QcTasks


workflow PlotVcfQcMetrics {
  input {
    Array[File] size_distribution_tsvs
    Array[File] af_distribution_tsvs
    Array[File] size_vs_af_distribution_tsvs
    File? all_svs_bed
    File? common_snvs_bed
    File? common_indels_bed
    File? common_svs_bed

    Float common_af_cutoff = 0.001

    String output_prefix

    String g2c_analysis_docker
  }

  ######################
  ### DATA PREPROCESSING
  ######################

  # If necessary, collapse size distribution into a single file
  if ( length(size_distribution_tsvs) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumSizeDistribs {
        input:
          distrib_tsvs = size_distribution_tsvs,
          out_prefix = output_prefix + ".size_distribution",
          g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File size_distrib = select_first(SumSizeDistribs.merged_distrib, size_distribution_tsvs[0]])

  # If necessary, collapse AF distribution into a single file
  if ( length(af_distribution_tsvs) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumAfDistribs {
        input:
          distrib_tsvs = af_distribution_tsvs,
          out_prefix = output_prefix + ".af_distribution",
          g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File af_distrib = select_first(SumAfDistribs.merged_distrib, af_distribution_tsvs[0]])

  # If necessary, collapse size-vs-AF distribution into a single file
  if ( length(size_vs_af_distribution_tsvs) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumJointDistribs {
        input:
          distrib_tsvs = size_vs_af_distribution_tsvs,
          out_prefix = output_prefix + ".size_vs_af_distribution",
          g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File joint_distrib = select_first(SumJointDistribs.merged_distrib, size_vs_af_distribution_tsvs[0]])

  #################
  ### VISUALIZATION
  #################

  # Generate site-level plots
  call PlotSiteMetrics {
    input:
      size_distrib = size_distrib,
      af_distrib = af_distrib,
      joint_distrib = joint_distrib,
      common_af_cutoff = common_af_cutoff,
      all_svs_bed = all_svs_bed,
      output_prefix = output_prefix,
      g2c_analysis_docker = g2c_analysis_docker
  }

  ##################
  ### OUTPUT CLEANUP
  ##################

  # Package all stats into a single tarball
  # TODO: implement this

  # Package all plots into a second, separate tarball
  # TODO: implement this

  output {}
}


task PlotSiteMetrics {
  input {
    File size_distrib
    File af_distrib
    File joint_distrib
    Float common_af_cutoff
    File? all_svs_bed
    File? common_snvs_bed
    File? common_indels_bed
    File? common_svs_bed

    String output_prefix

    Float mem_gb = 7.5
    Int n_cpu = 4
    Int? disk_gb

    String g2c_analysis_docker
  }

  Array[File?] loc_inputs = [size_distrib, af_distrib, joint_distrib, all_svs_bed,
                             common_snvs_bed, common_indels_bed, common_svs_bed]
  Int default_disk_gb = ceil(2 * size(select_all(loc_inputs), "GB")) + 20

  String pw_snv_cmd = if defined(common_snvs_bed) then "--snvs ~{common_snvs_bed}" else ""
  String pw_indel_cmd = if defined(common_indels_bed) then "--indels ~{common_indels_bed}" else ""
  String pw_sv_cmd = if defined(common_svs_bed) then "--indels ~{common_svs_bed}" else ""

  command <<<
    set -eu -o pipefail

    mkdir site_metrics

    # Plot site summary metrics
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_site_summary_metrics.R \
      --size-distrib ~{size_distrib} \
      --af-distrib ~{af_distrib} \
      --joint-distrib ~{joint_distrib} \
      --common-af ~{common_af_cutoff} \
      --out-prefix site_metrics/~{out_prefix}

    # Plot site-level metrics for common variants
    if ~{defined(common_snvs_bed)} || \
       ~{defined(common_indels_bed}) || \
       ~{defined(common_svs_bed)}; then
      /opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_site_pointwise_metrics.R \
        ~{pw_snv_cmd} \
        ~{pw_indel_cmd} \
        ~{pw_sv_cmd} \
        --combine \
        --out-prefix site_metrics/~{out_prefix}
    fi

    # Compress outputs
    tar -czvf site_metrics.tar.gz site_metrics/
  >>>

  output {
    File site_metric_plots_tarball = "site_metrics.tar.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 2
    max_retries: 1
  }

}
