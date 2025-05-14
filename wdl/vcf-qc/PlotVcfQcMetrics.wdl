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

    Array[File]? ref_size_distribution_tsvs
    Array[File]? ref_af_distribution_tsvs
    String ref_cohort_prefix = "ref_dataset"
    String ref_cohort_plot_title = "Ref. data"

    Array[File]? all_sv_beds
    Array[File]? common_snv_beds
    Array[File]? common_indel_beds
    Array[File]? common_sv_beds

    # Expected organization of site benchmarking inputs:
    # Outer array: one entry per benchmarking dataset
    # Middle arrays: one file per evaluation interval set
    # Inner arrays: one or more files, such as one file per chromosome 
    # (or just one file if already collapsed across all chromosomes)
    Array[Array[Array[File]]]? site_benchmark_common_snv_ppv_beds
    Array[Array[Array[File]]]? site_benchmark_common_indel_ppv_beds
    Array[Array[Array[File]]]? site_benchmark_common_sv_ppv_beds
    Array[Array[Array[File]]]? site_benchmark_ppv_by_freqs
    Array[Array[Array[File]]]? site_benchmark_sensitivity_by_freqs
    Array[String]? site_benchmark_dataset_prefixes
    Array[String]? site_benchmark_dataset_titles
    Array[String]? site_benchmark_interval_names

    Float common_af_cutoff = 0.001

    String output_prefix

    String bcftools_docker
    String g2c_analysis_docker
  }

  Boolean has_site_benchmarking = ( (defined(site_benchmark_common_snv_ppv_beds) 
                                     || defined(site_benchmark_common_indel_ppv_beds) 
                                     || defined(site_benchmark_common_sv_ppv_beds) 
                                     || defined(site_benchmark_ppv_by_freqs)
                                     || defined(site_benchmark_sensitivity_by_freqs)) 
                                    && defined(site_benchmark_dataset_prefixes) 
                                    && defined(site_benchmark_dataset_titles) 
                                    && defined(site_benchmark_interval_names))
  Int n_site_benchmark_datasets = length(select_first([site_benchmark_dataset_prefixes]))
  Int n_site_benchmark_intervals = length(select_first([site_benchmark_interval_names]))


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

  # If necessary, collapse reference size distribution into a single file
  if ( defined(ref_size_distribution_tsvs) ) {
    Array[File] rsdt_array = select_all(select_first([ref_size_distribution_tsvs]))
    if ( length(rsdt_array) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumRefSizeDistribs {
        input:
          distrib_tsvs = rsdt_array,
          out_prefix = ref_cohort_prefix + ".size_distribution",
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
    File ref_size_distrib = select_first([SumRefSizeDistribs.merged_distrib, rsdt_array[0]])
  }

  # If necessary, collapse reference AF distribution into a single file
  if ( defined(ref_af_distribution_tsvs) ) {
    Array[File] radt_array = select_all(select_first([ref_af_distribution_tsvs]))
    if ( length(radt_array) > 1 ) {
      call QcTasks.SumCompressedDistribs as SumRefAfDistribs {
        input:
          distrib_tsvs = radt_array,
          out_prefix = ref_cohort_prefix + ".af_distribution",
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
    File ref_af_distrib = select_first([SumRefAfDistribs.merged_distrib, radt_array[0]])
  }

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

  # Preprocess site benchmarking, if provided
  if (has_site_benchmarking) {
    scatter ( site_bench_di in range(n_site_benchmark_datasets) ) {

      String bd_name = flatten(select_all([site_benchmark_dataset_prefixes]))[site_bench_di]
      String sb_prefix = output_prefix + "." + bd_name

      scatter ( site_bench_ii in range(n_site_benchmark_intervals) ) {

        String bi_name = flatten(select_all([site_benchmark_interval_names]))[site_bench_ii]
        String sbi_prefix = sb_prefix + "." + bi_name
        
        # Collapse site benchmarking SNV BEDs
        if (defined(site_benchmark_common_snv_ppv_beds)) {
          call Utils.ConcatTextFiles as CollapseSiteBenchSnvs {
            input:
              shards = flatten(select_all([site_benchmark_common_snv_ppv_beds]))[site_bench_di][site_bench_ii],
              concat_command = "zcat",
              sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
              compression_command = "bgzip -c",
              input_has_header = true,
              output_filename = sbi_prefix + ".common_snvs.bed.gz",
              docker = bcftools_docker
          }
        }
        
        # Collapse site benchmarking indel BEDs
        if (defined(site_benchmark_common_indel_ppv_beds)) {
          call Utils.ConcatTextFiles as CollapseSiteBenchIndels {
            input:
              shards = flatten(select_all([site_benchmark_common_indel_ppv_beds]))[site_bench_di][site_bench_ii],
              concat_command = "zcat",
              sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
              compression_command = "bgzip -c",
              input_has_header = true,
              output_filename = sbi_prefix + ".common_indels.bed.gz",
              docker = bcftools_docker
          }
        }
        
        # Collapse site benchmarking SV BEDs
        if (defined(site_benchmark_common_sv_ppv_beds)) {
          call Utils.ConcatTextFiles as CollapseSiteBenchSvs {
            input:
              shards = flatten(select_all([site_benchmark_common_sv_ppv_beds]))[site_bench_di][site_bench_ii],
              concat_command = "zcat",
              sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
              compression_command = "bgzip -c",
              input_has_header = true,
              output_filename = sbi_prefix + ".common_svs.bed.gz",
              docker = bcftools_docker
          }
        }

        # Collapse compressed PPV tsvs
        if (defined(site_benchmark_ppv_by_freqs)) {
          call QcTasks.SumCompressedDistribs as SumSiteBenchPpvByAf {
            input:
              distrib_tsvs = select_first([site_benchmark_ppv_by_freqs])[site_bench_di][site_bench_ii],
              n_key_columns = 3,
              out_prefix = sbi_prefix + ".ppv_by_freq",
              g2c_analysis_docker = g2c_analysis_docker
          }
        }

        # Collapse compressed sensitivity tsvs
        if (defined(site_benchmark_sensitivity_by_freqs)) {
          call QcTasks.SumCompressedDistribs as SumSiteBenchSensByAf {
            input:
              distrib_tsvs = select_first([site_benchmark_sensitivity_by_freqs])[site_bench_di][site_bench_ii],
              n_key_columns = 3,
              out_prefix = sbi_prefix + ".sensitivity_by_freq",
              g2c_analysis_docker = g2c_analysis_docker
          }
        }
      }

      # For each benchmarking dataset, further collapse SNVs across interval sets
      if (defined(site_benchmark_common_snv_ppv_beds)) {
        call Utils.ConcatTextFiles as CollapseSiteBenchSnvsUnion {
            input:
              shards = select_all(CollapseSiteBenchSnvs.merged_file),
              concat_command = "zcat",
              sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
              compression_command = "bgzip -c",
              input_has_header = true,
              output_filename = sb_prefix + ".merged.common_snvs.bed.gz",
              docker = bcftools_docker
        }
      }

      # For each benchmarking dataset, further collapse indels across interval sets
      if (defined(site_benchmark_common_snv_ppv_beds)) {
        call Utils.ConcatTextFiles as CollapseSiteBenchIndelsUnion {
            input:
              shards = select_all(CollapseSiteBenchIndels.merged_file),
              concat_command = "zcat",
              sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
              compression_command = "bgzip -c",
              input_has_header = true,
              output_filename = sb_prefix + ".merged.common_indels.bed.gz",
              docker = bcftools_docker
        }
      }

      # For each benchmarking dataset, further collapse SVs across interval sets
      if (defined(site_benchmark_common_snv_ppv_beds)) {
        call Utils.ConcatTextFiles as CollapseSiteBenchSvsUnion {
            input:
              shards = select_all(CollapseSiteBenchSvs.merged_file),
              concat_command = "zcat",
              sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
              compression_command = "bgzip -c",
              input_has_header = true,
              output_filename = sb_prefix + ".merged.common_svs.bed.gz",
              docker = bcftools_docker
        }
      }
    }
  }
  Array[Array[File?]]? site_bench_snvs_merged = CollapseSiteBenchSnvs.merged_file
  Array[Array[File?]]? site_bench_indels_merged = CollapseSiteBenchIndels.merged_file
  Array[Array[File?]]? site_bench_svs_merged = CollapseSiteBenchSvs.merged_file
  Array[File]? site_bench_snvs_merged_union = select_all(select_first([CollapseSiteBenchSnvsUnion.merged_file]))
  Array[File]? site_bench_indels_merged_union = select_all(select_first([CollapseSiteBenchIndelsUnion.merged_file]))
  Array[File]? site_bench_svs_merged_union = select_all(select_first([CollapseSiteBenchSvsUnion.merged_file]))


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
      ref_size_distrib = ref_size_distrib,
      ref_af_distrib = ref_af_distrib,
      all_svs_bed = all_svs_bed,
      common_snvs_bed = common_snvs_bed,
      common_indels_bed = common_indels_bed,
      common_svs_bed = common_svs_bed,
      output_prefix = output_prefix,
      ref_title = ref_cohort_plot_title,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Plot site benchmarking, if provided
  if (has_site_benchmarking) {
    scatter ( site_bench_di in range(n_site_benchmark_datasets) ) {

      String plot_bd_prefix = select_first([site_benchmark_dataset_prefixes])[site_bench_di]
      String plot_bd_title = select_first([site_benchmark_dataset_titles])[site_bench_di]

      # Concatenate interval set-specific and union benchmarking BEDs for visualization
      Array[String] site_bench_di_names = flatten([select_first([site_benchmark_interval_names]), ["all"]])
      Array[File] site_bench_di_snv_to_plot = select_all(flatten([select_first([site_bench_snvs_merged])[site_bench_di], 
                                                                  [select_first([site_bench_snvs_merged_union])[site_bench_di]]]))
      Array[File] site_bench_di_indel_to_plot = select_all(flatten([select_first([site_bench_indels_merged])[site_bench_di], 
                                                                    [select_first([site_bench_indels_merged_union])[site_bench_di]]]))
      Array[File] site_bench_di_sv_to_plot = select_all(flatten([select_first([site_bench_svs_merged])[site_bench_di], 
                                                                 [select_first([site_bench_svs_merged_union])[site_bench_di]]]))
      Array[File] site_bench_di_ppv_by_af = select_all(select_first([SumSiteBenchPpvByAf.merged_distrib])[site_bench_di])
      Array[File] site_bench_di_sens_by_af = select_all(select_first([SumSiteBenchSensByAf.merged_distrib])[site_bench_di])

      call PlotSiteBenchmarking {
        input:
          ref_dataset_prefix = plot_bd_prefix,
          ref_dataset_title = plot_bd_title,
          eval_interval_names = site_bench_di_names,
          snv_beds = site_bench_di_snv_to_plot,
          indel_beds = site_bench_di_indel_to_plot,
          sv_beds = site_bench_di_sv_to_plot,
          ppv_by_af_tsvs = site_bench_di_ppv_by_af,
          sens_by_af_tsvs = site_bench_di_sens_by_af,
          output_prefix = output_prefix,
          common_af_cutoff = common_af_cutoff,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }


  ##################
  ### OUTPUT CLEANUP
  ##################

  # Package all stats into a single tarball
  # TODO: implement this

  # Package all plots into a second, separate tarball
  # TODO: implement this

  output {
    # For now, just outputting individual tarballs
    File site_metrics_tarball = PlotSiteMetrics.site_metric_plots_tarball
    # Array[File]? site_benchmarking_tarball = PlotSiteBenchmarking.site_benchmarking_plots_tarball
  }
}


task PlotSiteMetrics {
  input {
    File size_distrib
    File af_distrib
    File joint_distrib
    Float common_af_cutoff
    File? ref_size_distrib
    File? ref_af_distrib
    File? all_svs_bed
    File? common_snvs_bed
    File? common_indels_bed
    File? common_svs_bed

    String output_prefix
    String? ref_title

    Float mem_gb = 7.5
    Int n_cpu = 4
    Int? disk_gb

    String g2c_analysis_docker
  }

  Array[File?] loc_inputs = [size_distrib, af_distrib, joint_distrib, all_svs_bed,
                             common_snvs_bed, common_indels_bed, common_svs_bed]
  Int default_disk_gb = ceil(2 * size(select_all(loc_inputs), "GB")) + 20

  Boolean has_ref_size = defined(ref_size_distrib)
  String ref_size_bname = if has_ref_size then basename(select_first([ref_size_distrib])) else ""
  String ref_size_cmd = if has_ref_size then "--ref-size-distrib ~{ref_size_bname}" else ""

  Boolean has_ref_af = defined(ref_af_distrib)
  String ref_af_bname = if has_ref_af then basename(select_first([ref_af_distrib])) else ""
  String ref_af_cmd = if has_ref_af then "--ref-af-distrib ~{ref_af_bname}" else ""

  String ref_title_cmd = if defined(ref_title) then "--ref-title \"~{ref_title}\"" else ""

  Boolean has_all_svs = defined(all_svs_bed)
  String all_sv_bname = if has_all_svs then basename(select_first([all_svs_bed])) else ""
  String summary_sv_cmd = if has_all_svs then "--sv-sites ~{all_sv_bname}" else ""

  Boolean has_common_snvs = defined(common_snvs_bed)
  String common_snv_bname = if has_common_snvs then basename(select_first([common_snvs_bed])) else ""
  String pw_snv_cmd = if has_common_snvs then "--snvs ~{common_snv_bname}" else ""
  
  Boolean has_common_indels = defined(common_indels_bed)
  String common_indel_bname = if has_common_indels then basename(select_first([common_indels_bed])) else ""
  String pw_indel_cmd = if has_common_indels then "--indels ~{common_indel_bname}" else ""

  Boolean has_common_svs = defined(common_svs_bed)
  String common_sv_bname = if has_common_svs then basename(select_first([common_svs_bed])) else ""
  String pw_sv_cmd = if has_common_svs then "--svs ~{common_sv_bname}" else ""

  Boolean has_common_variants = has_common_snvs || has_common_indels || has_common_svs

  command <<<
    set -eu -o pipefail

    mkdir site_metrics

    # Symlink full SV BED to working directory
    if [ ~{has_all_svs} ]; then
      ln -s ~{default="" all_svs_bed} ~{all_sv_bname}
    fi

    # Symlink ref distribs to working directory
    if [ ~{has_ref_size} ]; then
      ln -s ~{default="" ref_size_distrib} ~{ref_size_bname}
    fi
    if [ ~{has_ref_af} ]; then
      ln -s ~{default="" ref_af_distrib} ~{ref_af_bname}
    fi

    # Plot site summary metrics
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_site_summary_metrics.R \
      --size-distrib ~{size_distrib} \
      --af-distrib ~{af_distrib} \
      --joint-distrib ~{joint_distrib} \
      ~{ref_size_cmd} \
      ~{ref_af_cmd} \
      ~{ref_title_cmd} \
      ~{summary_sv_cmd} \
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


task PlotSiteBenchmarking {
  input {
    String ref_dataset_prefix
    String ref_dataset_title
    Array[String] eval_interval_names

    Array[File]? snv_beds
    Array[File]? indel_beds
    Array[File]? sv_beds
    Array[File]? ppv_by_af_tsvs
    Array[File]? sens_by_af_tsvs

    String output_prefix

    Float common_af_cutoff

    Float mem_gb = 7.5
    Int n_cpu = 4
    Int? disk_gb

    String g2c_analysis_docker
  }

  Array[File] loc_inputs = flatten(select_all([snv_beds, indel_beds, sv_beds]))
  Int default_disk_gb = ceil(2 * size(loc_inputs, "GB")) + 20

  # Note that this outdir string is used as both a directory name and a file prefix
  String outdir = sub(output_prefix + "." + ref_dataset_prefix + "." + "site_benchmarking", "[ ]+", "_")

  Int n_sets = length(eval_interval_names)

  command <<<
    set -euo pipefail

    mkdir ~{outdir}

    # Write list of localized SNV beds
    if [ ~{defined(snv_beds)} ]; then
      cat ~{write_lines(select_first([snv_beds]))} > snv_beds.list
    else
      seq 1 ~{n_sets} | awk '{ print "." }' > snv_beds.list
    fi

    # Write list of localized indel beds
    if [ ~{defined(indel_beds)} ]; then
      cat ~{write_lines(select_first([indel_beds]))} > indel_beds.list
    else
      seq 1 ~{n_sets} | awk '{ print "." }' > indel_beds.list
    fi

    # Write list of localized SV beds
    if [ ~{defined(sv_beds)} ]; then
      cat ~{write_lines(select_first([sv_beds]))} > sv_beds.list
    else
      seq 1 ~{n_sets} | awk '{ print "." }' > sv_beds.list
    fi

    # Make input .tsv for pointwise plotting
    paste \
      ~{write_lines(eval_interval_names)} \
      snv_beds.list \
      indel_beds.list \
      sv_beds.list \
    > pointwise.in.tsv
    echo "Pointwise benchmarking inputs parsed as follows:"
    cat pointwise.in.tsv

    # Perform pointwise plotting once for each eval interval set
    while read set_name snv_bed indel_bed sv_bed; do
      cmd="/opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_pointwise_site_benchmarking.R"
      if [ $snv_bed != "." ]; then
        cmd="$cmd --snvs $snv_bed"
      fi
      if [ $indel_bed != "." ]; then
        cmd="$cmd --indels $indel_bed"
      fi
      if [ $sv_bed != "." ]; then
        cmd="$cmd --svs $sv_bed"
      fi
      cmd="$cmd --common-af ~{common_af_cutoff} --ref-title \"~{ref_dataset_title}\""
      cmd="$cmd --combine --out-prefix ~{outdir}/~{outdir}.$set_name"
      echo -e "Now performing pointwise site benchmarking as follows:\n$cmd"
      eval "$cmd"
    done < pointwise.in.tsv

    # Plot site summary metrics like PPV and sensitivity
    cmd="/opt/pancan_germline_wgs/scripts/qc/vcf_qc/plot_site_benchmarking_metrics.R"
    if [ ~{defined(ppv_by_af_tsvs)} == "true" ]; then
      while read tsv; do
        cmd="$cmd --ppv-by-af $tsv"
      done < ~{write_lines(select_first([ppv_by_af_tsvs]))}
    fi
    if [ ~{defined(sens_by_af_tsvs)} == "true" ]; then
      while read tsv; do
        cmd="$cmd --sens-by-af $tsv"
      done < ~{write_lines(select_first([sens_by_af_tsvs]))}
    fi
    while read sname; do
      cmd="$cmd --set-name $sname"
    done < ~{write_lines(eval_interval_names)}
    cmd="$cmd --ref-title \"~{ref_dataset_title}\" --common-af ~{common_af_cutoff}"
    cmd="$cmd --out-prefix ~{outdir}/~{outdir}"
    echo -e "Now performing site benchmarking metric visualization as follows:\n$cmd"
    eval "$cmd"

    # Compress outputs
    tar -czvf "~{outdir}.tar.gz" \
      ~{outdir}/
  >>>

  output {
    File site_benchmarking_plots_tarball = "~{outdir}.tar.gz"
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

