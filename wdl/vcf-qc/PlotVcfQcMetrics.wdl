# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot quality metrics for assessing a GATK joint-genotyped VCF (SNVs/indels or SVs)
# See CollectVcfQcMetrics.wdl for the generation of inputs for this workflow


version 1.0


import "QcTasks.wdl" as QcTasks
import "Utilities.wdl" as Utils


workflow PlotVcfQcMetrics {
  input {
    Array[File] size_distribution_tsvs
    Array[File] af_distribution_tsvs
    Array[File] size_vs_af_distribution_tsvs

    Array[File]? all_sv_beds
    Array[File]? common_snv_beds
    Array[File]? common_indel_beds
    Array[File]? common_sv_beds

    Float common_af_cutoff = 0.001

    String output_prefix

    String bcftools_docker
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
  File size_distrib = select_first([SumSizeDistribs.merged_distrib, size_distribution_tsvs[0]])

  # If necessary, collapse AF distribution into a single file
  if ( length(af_distribution_tsvs) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumAfDistribs {
        input:
          distrib_tsvs = af_distribution_tsvs,
          out_prefix = output_prefix + ".af_distribution",
          g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File af_distrib = select_first([SumAfDistribs.merged_distrib, af_distribution_tsvs[0]])

  # If necessary, collapse size-vs-AF distribution into a single file
  if ( length(size_vs_af_distribution_tsvs) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumJointDistribs {
        input:
          distrib_tsvs = size_vs_af_distribution_tsvs,
          n_key_columns = 3,
          out_prefix = output_prefix + ".size_vs_af_distribution",
          g2c_analysis_docker = g2c_analysis_docker
    }
  }
  File joint_distrib = select_first([SumJointDistribs.merged_distrib, size_vs_af_distribution_tsvs[0]])

  # If necessary, collapse all SV BEDs
  if ( length(select_first(select_all([all_sv_beds]))) > 1 ) {
    call Utils.ConcatTextFiles as CollapseAllSvs {
      input:
        shards = select_first([all_sv_beds]),
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".all_svs.bed.gz",
        docker = bcftools_docker
      }
  }
  File? all_svs_bed = select_first(select_all([CollapseAllSvs.merged_file, 
                                               select_first([all_sv_beds])[0]]))

  # If necessary, collapse common SNV BEDs
  if ( length(select_first(select_all([common_snv_beds]))) > 1 ) {
    call Utils.ConcatTextFiles as CollapseCommonSnvs {
      input:
        shards = select_first([common_snv_beds]),
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".common_snvs.bed.gz",
        docker = bcftools_docker
      }
  }
  File? common_snvs_bed = select_first(select_all([CollapseCommonSnvs.merged_file, 
                                                   select_first([common_snv_beds])[0]]))

  # If necessary, collapse common indel BEDs
  if ( length(select_first(select_all([common_indel_beds]))) > 1 ) {
    call Utils.ConcatTextFiles as CollapseCommonIndels {
      input:
        shards = select_first([common_indel_beds]),
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".common_indels.bed.gz",
        docker = bcftools_docker
      }
  }
  File? common_indels_bed = select_first(select_all([CollapseCommonIndels.merged_file, 
                                                     select_first([common_indel_beds])[0]]))

  # If necessary, collapse common SV BEDs
  if ( length(select_first(select_all([common_sv_beds]))) > 1 ) {
    call Utils.ConcatTextFiles as CollapseCommonSvs {
      input:
        shards = select_first([common_sv_beds]),
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".common_svs.bed.gz",
        docker = bcftools_docker
      }
  }
  File? common_svs_bed = select_first(select_all([CollapseCommonSvs.merged_file, 
                                                  select_first([common_sv_beds])[0]]))

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
      common_snvs_bed = common_snvs_bed,
      common_indels_bed = common_indels_bed,
      common_svs_bed = common_svs_bed,
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

  output {
    # For now, just outputting site metrics tarball
    File site_metrics_tarball = PlotSiteMetrics.site_metric_plots_tarball
  }
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

  Boolean has_common_snvs = defined(common_snvs_bed)
  String common_snv_bname = if has_common_snvs then basename(select_first([common_snvs_bed])) else ""
  String pw_snv_cmd = if has_common_snvs then "--snvs common_snv_bname" else ""
  
  Boolean has_common_indels = defined(common_indels_bed)
  String common_indel_bname = if has_common_indels then basename(select_first([common_indels_bed])) else ""
  String pw_indel_cmd = if has_common_indels then "--indels common_indel_bname" else ""

  Boolean has_common_svs = defined(common_svs_bed)
  String common_sv_bname = if has_common_svs then basename(select_first([common_svs_bed])) else ""
  String pw_sv_cmd = if has_common_svs then "--svs common_sv_bname" else ""

  Boolean has_common_variants = has_common_snvs || has_common_indels || has_common_svs

  command <<<
    set -eu -o pipefail

    mkdir site_metrics

    # Plot site summary metrics
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_site_summary_metrics.R \
      --size-distrib ~{size_distrib} \
      --af-distrib ~{af_distrib} \
      --joint-distrib ~{joint_distrib} \
      --common-af ~{common_af_cutoff} \
      --out-prefix site_metrics/~{output_prefix}

    # Symlink common variant BEDs to working directory
    if [ ~{has_common_snvs} ]; then
      ln -s ~{default="" common_snvs_bed} ~{common_snv_bname}
    fi
    if [ ~{has_common_indels} ]; then
      ln -s ~{default="" common_indels_bed} ~{common_indel_bname}
    fi
    if [ ~{has_common_svs} ]; then
      ln -s ~{default="" common_svs_bed} ~{common_sv_bname}
    fi

    # Plot site-level metrics for common variants
    if ~{has_common_variants}; then
      /opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_site_pointwise_metrics.R \
        ~{pw_snv_cmd} \
        ~{pw_indel_cmd} \
        ~{pw_sv_cmd} \
        --combine \
        --out-prefix site_metrics/~{output_prefix}
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
