# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collect quality metrics for assessing a GATK joint-genotyped VCF (SNVs/indels and/or SVs)


version 1.0


import "BenchmarkSites.wdl" as BenchSites
import "QcTasks.wdl" as QcTasks
import "Utilities.wdl" as Utils


workflow CollectVcfQcMetrics {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    File? trios_fam_file                        # .fam file of known families for Mendelian transmission analyses
    File? sample_priority_list                  # Rank-ordered list of samples to retain for sample-level analyses
    Int n_for_sample_level_analyses = 1000      # Number of samples to use for sample-level summary analyses

    Boolean shard_vcf = true                    # Should the input VCF be sharded for QC collection?
    File? scatter_intervals_list                # GATK-style intervals file for scattering over vcf 
                                                # (any tabix-compliant interval definitions should work)
    Int n_records_per_shard = 25000             # Number of records per shard. This will only be used as a backup if 
                                                # scatter_intervals_list is not provided and shard_vcf is true

    Float common_af_cutoff = 0.001              # Minimum AF for a variant to be included in common variant subsets

    Array[File?] snv_site_benchmark_beds        # BED files for SNV site benchmarking; one per reference dataset or cohort
    Array[File?] indel_site_benchmark_beds      # BED files for SNV site benchmarking; one per reference dataset or cohort
    Array[File?] sv_site_benchmark_beds         # BED files for SNV site benchmarking; one per reference dataset or cohort
    Array[String?] site_benchmark_dataset_names

    Array[File]? benchmark_interval_beds        # BED files of intervals to consider for benchmarking evaluation
    Array[String]? benchmark_interval_bed_names # Descriptive names for each set of evaluation intervals
    File? genome_file                           # BEDTools-style .genome file
    Int benchmarking_shards = 2500              # Number of parallel tasks to use for site benchmarking

    String output_prefix

    String bcftools_docker
    String g2c_analysis_docker
    String linux_docker
  }

  #####################
  ### SAMPLE MANAGEMENT
  #####################

  # Get list of samples present in input VCFs
  # Uses the last input VCF, assuming it will be the smallest in most cases
  Int last_vcf_index = length(vcfs) - 1
  call Utils.GetSamplesFromVcfHeader as GetSamplesInVcf {
    input:
      vcf = vcfs[last_vcf_index],
      vcf_idx = vcf_idxs[last_vcf_index],
      bcftools_docker = bcftools_docker
  }

  # Clean input .fam file to retain only complete trios
  if( defined(trios_fam_file) ) {
    call CleanFam {
      input:
        fam_file = select_first([trios_fam_file]),
        all_samples_list = GetSamplesInVcf.sample_list,
        docker = g2c_analysis_docker
    }
  }
  File all_qc_samples_list = select_first([CleanFam.all_samples_no_probands_list,
                                           GetSamplesInVcf.sample_list])

  # Define list of samples to use for sample-specific analyses
  call ChooseTargetSamples {
    input:
      all_samples_list = all_qc_samples_list,
      sample_priority_list = sample_priority_list,
      n_samples = n_for_sample_level_analyses,
      out_prefix = output_prefix,
      docker = bcftools_docker
  }

  #####################
  ### VCF PREPROCESSING
  #####################

  # Read scatter intervals, if optioned
  if ( defined(scatter_intervals_list) ) {
    call Utils.ParseIntervals {
      input:
        intervals_list = select_first([scatter_intervals_list]),
        docker = linux_docker
    }
  }

  # Preprocess each VCF according to desired scatter behavior
  Array[Pair[File, File]] input_vcf_infos = zip(vcfs, vcf_idxs)
  scatter ( input_vcf_info in input_vcf_infos ) {

    File vcf = input_vcf_info.left
    File vcf_idx = input_vcf_info.right

    # Check the header of each input VCF for the presence of mCNVs
    call Utils.McnvHeaderCheck as McnvCheck {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        docker = g2c_analysis_docker
    }

    if ( shard_vcf ) {

      # Default to scattering over prespecified intervals, as this is most compute efficient
      if ( defined(scatter_intervals_list) ) {        
        scatter ( interval_info in select_first([ParseIntervals.interval_info]) ) {

          String interval_coords = interval_info[1]
          String shard_prefix = basename(vcf, ".vcf.gz") + "." + interval_info[0]

          # Extract desired interval from VCF
          call Utils.StreamSliceVcf as SliceInterval {
            input:
              vcf = vcf,
              vcf_idx = vcf_idx, 
              interval = interval_coords,
              outfile_name = shard_prefix + ".vcf.gz",
              bcftools_docker = bcftools_docker
          }
        }
      }

      # If scatter_intervals_list is not provided, shard the VCF the old fashioned way
      if ( !defined(scatter_intervals_list) ) {

        # Shard VCF on a VM (costly & slow for large VCFs)
        call Utils.ShardVcf {
          input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            records_per_shard = n_records_per_shard,
            bcftools_docker = bcftools_docker,
            n_preemptible = 1
        }
      }
    }

    # Unify the three approaches for sharding a VCF prior to preprocessing
    Array[File] vcf_shards = select_first([SliceInterval.vcf_slice,
                                           ShardVcf.vcf_shards,
                                           [vcf]])
    Array[File] vcf_shard_idxs = select_first([SliceInterval.vcf_slice_idx,
                                               ShardVcf.vcf_shard_idxs,
                                               [vcf_idx]])

    # Scatter over VCF shards and preprocess each shard
    Array[Pair[File, File]] pp_vcf_infos = zip(vcf_shards, vcf_shard_idxs)
    scatter ( pp_vcf_info in pp_vcf_infos ) {

      File pp_vcf = pp_vcf_info.left
      File pp_vcf_idx = pp_vcf_info.right

      call PreprocessVcf {
        input:
          vcf = pp_vcf,
          vcf_idx = pp_vcf_idx,
          target_samples = ChooseTargetSamples.target_samples,
          trio_samples = CleanFam.trio_samples_list,
          has_mcnvs = McnvCheck.has_mcnvs,
          out_prefix = basename(vcf, ".vcf.gz"),
          docker = bcftools_docker
      }
    }
  }

  ###########################
  ### BASIC METRIC COLLECTION
  ###########################

  Array[File] site_vcf_shards = flatten(PreprocessVcf.sites_vcf)
  Array[File] site_vcf_shard_idxs = flatten(PreprocessVcf.sites_vcf_idx)
  Array[Pair[File, File]] site_vcf_info = zip(site_vcf_shards, site_vcf_shard_idxs)

  # Collect site-level metrics for all preprocessed shards
  scatter ( shard_info in site_vcf_info ) {

    # Compute site-level metrics
    call QcTasks.CollectSiteMetrics {
      input:
        vcf = shard_info.left,
        vcf_idx = shard_info.right,
        n_samples = ChooseTargetSamples.n_samples_all,
        common_af_cutoff = common_af_cutoff,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Compute sample-level metrics
    # TODO: implement this

    # Compute trio metrics
    # TODO: implement this

    # Compute sample-level benchmarking metrics
    # TODO: implement this
  }

  # Collapse all SNV sites
  Array[File] all_snv_beds = select_all(CollectSiteMetrics.snv_sites)
  if ( length(all_snv_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseAllSnvs {
      input:
        shards = all_snv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".all_snvs.bed.gz",
        docker = bcftools_docker
      }
  }

  # Collapse common SNV sites
  Array[File] common_snv_beds = select_all(CollectSiteMetrics.common_snv_sites)
  if ( length(common_snv_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseCommonSnvs {
      input:
        shards = common_snv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".common_snvs.bed.gz",
        docker = bcftools_docker
      }
  }

  # Collapse all indel sites
  Array[File] all_indel_beds = select_all(CollectSiteMetrics.indel_sites)
  if ( length(all_indel_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseAllIndels {
      input:
        shards = all_indel_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".all_indels.bed.gz",
        docker = bcftools_docker
      }
  }

  # Collapse common indel sites
  Array[File] common_indel_beds = select_all(CollectSiteMetrics.common_indel_sites)
  if ( length(common_indel_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseCommonIndels {
      input:
        shards = common_indel_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".common_indels.bed.gz",
        docker = bcftools_docker
      }
  }

  # Collapse all SV sites
  Array[File] all_sv_beds = select_all(CollectSiteMetrics.sv_sites)
  if ( length(all_sv_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseAllSvs {
      input:
        shards = all_sv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".all_svs.bed.gz",
        docker = bcftools_docker
      }
  }

  # Collapse common SV sites
  Array[File] common_sv_beds = select_all(CollectSiteMetrics.common_sv_sites)
  if ( length(common_sv_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseCommonSvs {
      input:
        shards = common_sv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".common_svs.bed.gz",
        docker = bcftools_docker
      }
  }

  #########################
  ### EXTERNAL BENCHMARKING
  #########################

  # Compute site-level benchmarking metrics
  Int n_site_benchmark_datasets = length(site_benchmark_dataset_names)

  if ( n_site_benchmark_datasets > 0 ) {
    scatter ( site_bench_idx in range(n_site_benchmark_datasets) ) {
      call BenchSites.BenchmarkSites {
        input:
          source_snv_bed = CollapseAllSnvs.merged_file,
          source_indel_bed = CollapseAllIndels.merged_file,
          source_sv_bed = CollapseAllSvs.merged_file,
          source_prefix = output_prefix,
          target_snv_bed = select_all(snv_site_benchmark_beds)[site_bench_idx],
          target_indel_bed = select_all(indel_site_benchmark_beds)[site_bench_idx],
          target_sv_bed = select_all(sv_site_benchmark_beds)[site_bench_idx],
          target_prefix = select_all(site_benchmark_dataset_names)[site_bench_idx],
          eval_interval_beds = select_first([benchmark_interval_beds]),
          eval_interval_bed_names = select_first([benchmark_interval_bed_names]),
          genome_file = select_first([genome_file]),
          total_shards = benchmarking_shards,
          common_af_cutoff = common_af_cutoff,
          bcftools_docker = bcftools_docker,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }
  

  ##################
  ### OUTPUT CLEANUP
  ##################

  # Collapse size distributions
  call QcTasks.SumCompressedDistribs as SumSizeDistribs {
    input:
      distrib_tsvs = CollectSiteMetrics.size_distrib,
      out_prefix = output_prefix + ".size_distribution",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collapse freq distributions
  call QcTasks.SumCompressedDistribs as SumAfDistribs {
    input:
      distrib_tsvs = CollectSiteMetrics.af_distrib,
      out_prefix = output_prefix + ".af_distribution",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collapse size X freq distributions
  call QcTasks.SumCompressedDistribs as SumJointDistribs {
    input:
      distrib_tsvs = CollectSiteMetrics.size_vs_af_distrib,
      n_key_columns = 3,
      out_prefix = output_prefix + ".size_vs_af_distribution",
      g2c_analysis_docker = g2c_analysis_docker
  }

  output {
    File? all_snvs_bed = CollapseAllSnvs.merged_file
    File? all_indels_bed = CollapseAllIndels.merged_file
    File? all_svs_bed = CollapseAllSvs.merged_file
    File? common_snvs_bed = CollapseCommonSnvs.merged_file
    File? common_indels_bed = CollapseCommonIndels.merged_file
    File? common_svs_bed = CollapseCommonSvs.merged_file

    File size_distrib = SumSizeDistribs.merged_distrib
    File af_distrib = SumAfDistribs.merged_distrib
    File size_vs_af_distrib = SumJointDistribs.merged_distrib

    Array[Array[File?]]? site_benchmark_common_snv_ppv_beds = BenchmarkSites.common_snv_ppv_beds
    Array[Array[File?]]? site_benchmark_common_snv_sens_beds = BenchmarkSites.common_snv_sens_beds
    Array[Array[File?]]? site_benchmark_common_indel_ppv_beds = BenchmarkSites.common_indel_ppv_beds
    Array[Array[File?]]? site_benchmark_common_indel_sens_beds = BenchmarkSites.common_indel_sens_beds
    Array[Array[File?]]? site_benchmark_common_sv_ppv_beds = BenchmarkSites.common_sv_ppv_beds
    Array[Array[File?]]? site_benchmark_common_sv_sens_beds = BenchmarkSites.common_sv_sens_beds
    Array[Array[File]]? site_benchmark_ppv_by_sizes = BenchmarkSites.ppv_by_sizes
    Array[Array[File]]? site_benchmark_sensitivity_by_sizes = BenchmarkSites.sensitivity_by_sizes
    Array[Array[File]]? site_benchmark_ppv_by_freqs = BenchmarkSites.ppv_by_freqs
    Array[Array[File]]? site_benchmark_sensitivity_by_freqs = BenchmarkSites.sensitivity_by_freqs
  }
}


# Clean input .fam file to restrict to complete trios present in input VCF
task CleanFam {
  input {
    File fam_file
    File all_samples_list
    String docker
  }

  String fam_out_fname = basename(fam_file, ".fam") + ".cleaned_trios.fam"
  String proband_out_fname = basename(fam_file, ".fam") + ".proband_ids.list"
  String unrelated_out_fname = basename(all_samples_list) + ".no_probands.list"
  String trio_samples_out_fname = basename(fam_file, ".fam") + ".trio_samples.list"
  
  command <<<
    set -eu -o pipefail

    # Clean .fam file to duos or trios present in all_samples_list
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/subset_fam.R \
      --in-fam ~{fam_file} \
      --all-samples ~{all_samples_list} \
      --out-fam duos_and_trios.fam

    # Write list of all probands (can exclude these samples for better AF estimates)
    cut -f1 duos_and_trios.fam > "~{proband_out_fname}"

    # Remove probands from list of all samples (same rationale as above)
    fgrep \
      -xvf "~{proband_out_fname}" \
      ~{all_samples_list} \
    > "~{unrelated_out_fname}"

    # Further restrict .fam to complete trios
    awk -v FS="\t" -v OFS="\t" \
      '{ if ($3!=0 && $4!=0) print }' \
      duos_and_trios.fam \
    > "~{fam_out_fname}"

    # Write list of all samples in complete trios
    awk -v OFS="\n" '{ print $2, $3, $4 }' "~{fam_out_fname}" \
    | sort -V > "~{trio_samples_out_fname}"
  >>>

  output {
    File trios_fam = "~{fam_out_fname}"
    File probands_list = "~{proband_out_fname}"
    File all_samples_no_probands_list = "~{unrelated_out_fname}"
    File trio_samples_list = "~{trio_samples_out_fname}"
  }

  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk 20 HDD"
    preemptible: 3
  }
}


# Define list of "target" samples to use for sample-based analyses
task ChooseTargetSamples {
  input {
    File all_samples_list
    String out_prefix
    File? sample_priority_list
    Int n_samples
    String docker
  }

  String outfile = out_prefix + ".target.samples.list"

  command <<<
    set -eu -o pipefail

    touch "~{outfile}"

    # First, add samples from a priority list, if optioned
    if [ ~{defined(sample_priority_list)} == "true" ]; then
      fgrep \
        -xf ~{all_samples_list} \
        ~{select_first([sample_priority_list, ""])} \
      | head -n ~{n_samples} \
      >> "~{outfile}" || true
      fgrep \
        -xvf ~{select_first([sample_priority_list, ""])} \
        ~{all_samples_list} \
      > remainder.samples.list
    else
      cp ~{all_samples_list} remainder.samples.list
    fi

    # Second, supplement from random selection of remaining samples
    subtotal=$( cat ~{outfile} | wc -l )
    n_to_add=$(( ~{n_samples} - $subtotal ))
    if [ $n_to_add -gt 0 ]; then
      shuf \
        --random-source=<( yes 2025 ) \
        remainder.samples.list \
      | head -n $n_to_add \
      >> "~{outfile}" || true
    fi
  >>>

  output {
    File target_samples = outfile
    Int n_samples_all = length(read_lines(all_samples_list))
  }

  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk 20 HDD"
    preemptible: 3
  }
}


# Update all VCF INFO required for QC metric collection, make all necessary genotype subsets
task PreprocessVcf {
  input {
    File vcf
    File vcf_idx
    File target_samples
    File? trio_samples
    File? benchmarking_samples
    Boolean has_mcnvs = false
    
    String out_prefix
    
    Int? disk_gb
    Float mem_gb = 4
    Int n_cpu = 2
    
    String docker
  }

  Array[File] all_sample_lists = select_all([trio_samples, target_samples, benchmarking_samples])

  String sites_outfile = out_prefix + ".sites.vcf.gz"

  String mcnv_anno = if has_mcnvs then "| /opt/pancan_germline_wgs/scripts/gatksv_helpers/annotate_mcnv_freqs.py - -" else ""

  Int default_disk_gb = ceil(3 * size(vcf, "GB")) + 10
  Int hdd_gb = select_first([disk_gb, default_disk_gb])

  command <<<
    set -eu -o pipefail

    # Define superset of samples whose GTs we need for subsequent analysis
    cat ~{sep=" " all_sample_lists} | sort -V | uniq > all.samples.list

    # Preprocess VCF
    bcftools +fill-tags ~{vcf} -- -t AN,AC,AF,AC_Hemi,AC_Het,AC_Hom,HWE \
    ~{mcnv_anno} \
    | bcftools view \
      --samples-file all.samples.list \
      --no-update \
    -Oz -o "~{out_prefix}.cleaned_wGTs.vcf.gz"
    tabix -p vcf "~{out_prefix}.cleaned_wGTs.vcf.gz"

    # Make sites VCF for site-level analyses
    bcftools view -G \
      -Oz -o "~{sites_outfile}" \
      "~{out_prefix}.cleaned_wGTs.vcf.gz"
    tabix -p vcf "~{sites_outfile}"

    # Make VCF of only sites & samples appearing in target samples
    # TODO: implement this

    # Make VCF of only sites & samples appearing in complete trios
    # TODO: implement this

    # Make VCF of only sites & samples appearing in benchmarking samples
    # TODO: implement this
  >>>

  output {
    File sites_vcf = sites_outfile
    File sites_vcf_idx = "~{sites_outfile}.tbi"
  }

  runtime {
    docker: docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk ~{hdd_gb} HDD"
    preemptible: 3
  }
}

