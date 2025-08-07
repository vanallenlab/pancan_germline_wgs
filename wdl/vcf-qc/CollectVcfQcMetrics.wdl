# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collect quality metrics for assessing a GATK joint-genotyped VCF (SNVs/indels and/or SVs)


version 1.0


import "BenchmarkSamples.wdl" as BenchSamples
import "BenchmarkSites.wdl" as BenchSites
import "BenchmarkTrios.wdl" as BenchTrios
import "BenchmarkTwins.wdl" as BenchTwins
import "QcTasks.wdl" as QcTasks


workflow CollectVcfQcMetrics {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    Boolean shard_vcf = true                       # Should the input VCF be sharded for QC collection?
    File? scatter_intervals_list                   # GATK-style intervals file for scattering over vcf 
                                                   # (any tabix-compliant interval definitions should work)
    Int n_records_per_shard = 25000                # Number of records per shard. This will only be used as a backup if 
                                                   # scatter_intervals_list is not provided and shard_vcf is true
    String extra_vcf_preprocessing_commands = ""   # Optional string of extra shell commands to execute when preprocessing
                                                   # main input vcfs. This must be prefixed by a pipe so that it can
                                                   # be injected into a chain of bash commands while reading from stdin 
                                                   # and writing to stdout

    Float common_af_cutoff = 0.01                  # Minimum AF for a variant to be included in common variant subsets

    File? trios_fam_file                           # .fam file of trios for Mendelian transmission analyses
    File? twins_tsv                                # Two-column .tsv with pairs of IDs for technical replicates or identical twins
    File? sample_priority_tsv                      # Three-column .tsv of sample ID, priority tier (integer), and sampling weight (float)
                                                   # This file will be used in conjunction with n_for_sample_level_analyses
                                                   # to determine which samples will be retained for sample-level analyses. Note that
                                                   # this will also influence the number of samples included for sample-level 
                                                   # benchmarking and twin/trio analysis. If you want to enforce a certain set of 
                                                   # samples to always be included for these analyses, you should specify their IDs
                                                   # with high priority tier in this file; otherwise, sample overlaps will be 
                                                   # left to random chance, which will often be suboptimal.
    Int n_for_sample_level_analyses = 1000         # Number of samples to use for all sample-level analyses, including trio/twin/benchmarking

    Array[File?] snv_site_benchmark_beds           # BED files for SNV site benchmarking; one per reference dataset or cohort
    Array[File?] indel_site_benchmark_beds         # BED files for SNV site benchmarking; one per r`eference dataset or cohort
    Array[File?] sv_site_benchmark_beds            # BED files for SNV site benchmarking; one per reference dataset or cohort
    Array[String?] site_benchmark_dataset_names

    Array[Array[File?]] sample_benchmark_vcfs      # VCFs to use for sample-level genotype benchmarking. Each outer array corresponds
                                                   # to a single dataset or technology/modality. Each inner array can have any number
                                                   # of VCFs; each of these VCFs will be benchmarked in parallel and their results
                                                   # will be concatenated within each sample
    Array[Array[File?]] sample_benchmark_vcf_idxs  # Tabix indexes for each VCF in sample_benchmark_vcfs; order must match.
    Array[Array[File?]] sample_benchmark_id_maps   # Two-column .tsv mapping IDs in main input VCF to each VCF in sample_benchmark_vcfs
                                                   # One id_map must be provided for each VCF in sample_benchmark_vcfs
    Array[String?] sample_benchmark_dataset_names

    Array[File]? benchmark_interval_beds           # BED files of intervals to consider for benchmarking evaluation
    Array[String]? benchmark_interval_bed_names    # Descriptive names for each set of evaluation intervals
    File? genome_file                              # BEDTools-style .genome file
    Int benchmarking_shards = 2500                 # Number of parallel tasks to use for site and sample benchmarking
    Int min_samples_per_bench_shard = 10           # Minimum number of samples per shard to allow for sample benchmarking

    String output_prefix

    String bcftools_docker
    String g2c_analysis_docker
    String linux_docker
  }

  Int n_site_benchmark_datasets = length(site_benchmark_dataset_names)
  Boolean do_site_bench = (n_site_benchmark_datasets > 0)

  Int n_gt_benchmark_datasets = length(sample_benchmark_dataset_names)
  Boolean do_sample_bench = (n_gt_benchmark_datasets > 0)

  #####################
  ### SAMPLE MANAGEMENT
  #####################

  # Get list of samples present in input VCFs
  # Uses the last input VCF, assuming it will be the smallest in most cases
  Int last_vcf_index = length(vcfs) - 1
  call QcTasks.GetSamplesFromVcfHeader as GetSamplesInVcf {
    input:
      vcf = vcfs[last_vcf_index],
      vcf_idx = vcf_idxs[last_vcf_index],
      bcftools_docker = bcftools_docker
  }

  # Clean input .fam file to retain only complete trios
  if ( defined(trios_fam_file) ) {
    call CleanFam {
      input:
        fam_file = select_first([trios_fam_file]),
        all_samples_list = GetSamplesInVcf.sample_list,
        docker = g2c_analysis_docker
    }
  }
  Boolean has_trios = select_first([CleanFam.n_trios, 0]) > 0

  # Clean input twins file to retain only complete pairs
  if ( defined(twins_tsv) ) {
    call CleanTwins {
      input:
        twins_tsv = select_first([twins_tsv]),
        all_samples_list = GetSamplesInVcf.sample_list,
        docker = g2c_analysis_docker
    }
  }
  Boolean has_twins = select_first([CleanTwins.n_twins, 0]) > 0

  # Define list of samples to use for sample-specific analyses
  Int actual_n_for_sample_level = if n_for_sample_level_analyses > GetSamplesInVcf.n_samples 
                                  then GetSamplesInVcf.n_samples 
                                  else n_for_sample_level_analyses
  call ChooseTargetSamples {
    input:
      all_samples_list = GetSamplesInVcf.sample_list,
      probands_list = CleanFam.probands_list,
      duplicate_samples_list = CleanTwins.duplicate_samples_list,
      sample_priority_tsv = sample_priority_tsv,
      n_samples = actual_n_for_sample_level,
      out_prefix = output_prefix,
      g2c_analysis_docker = g2c_analysis_docker
  }

  #####################
  ### VCF PREPROCESSING
  #####################

  call QcTasks.MakeHeaderFiller {}

  # Read scatter intervals, if optioned
  if ( defined(scatter_intervals_list) ) {
    call QcTasks.ParseIntervals {
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
    call QcTasks.McnvHeaderCheck as McnvCheck {
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
          call QcTasks.StreamSliceVcf as SliceInterval {
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
        call QcTasks.ShardVcf {
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
          site_exclude_samples = ChooseTargetSamples.site_exclude_samples,
          has_mcnvs = McnvCheck.has_mcnvs,
          extra_commands = extra_vcf_preprocessing_commands,
          supp_vcf_header = MakeHeaderFiller.supp_vcf_header,
          out_prefix = basename(vcf, ".vcf.gz"),
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Declare or create dummy files for default fills later  
  call QcTasks.MakeDummyFile as MakeEmptyFile {}
  call QcTasks.MakeEmptyBenchBed {
    input:
      docker = bcftools_docker
  }
  File empty_vcf = flatten(PreprocessVcf.empty_vcf)[0]
  File empty_vcf_idx = flatten(PreprocessVcf.empty_vcf_idx)[0]
  File empty_bed = MakeEmptyBenchBed.empty_bed

  ###########################
  ### BASIC METRIC COLLECTION
  ###########################

  Array[File] site_vcf_shards = flatten(PreprocessVcf.sites_vcf)
  Array[File] site_vcf_shard_idxs = flatten(PreprocessVcf.sites_vcf_idx)
  Array[Pair[File, File]] site_vcf_info = zip(site_vcf_shards, site_vcf_shard_idxs)

  Array[File] dense_vcf_shards = flatten(PreprocessVcf.dense_vcf)
  Array[File] dense_vcf_shard_idxs = flatten(PreprocessVcf.dense_vcf_idx)
  Array[Pair[File, File]] dense_vcf_info = zip(dense_vcf_shards, dense_vcf_shard_idxs)

  Array[File] dense_sites_vcf_shards = flatten(PreprocessVcf.dense_sites_vcf)
  Array[File] dense_sites_vcf_shard_idxs = flatten(PreprocessVcf.dense_sites_vcf_idx)
  Array[Pair[File, File]] dense_sites_vcf_info = zip(dense_sites_vcf_shards, dense_sites_vcf_shard_idxs)

  # Collect site-level metrics for all preprocessed shards
  scatter ( i in range(length(site_vcf_info)) ) {

    File sites_vcf_shard = site_vcf_info[i].left
    File sites_vcf_shard_idx = site_vcf_info[i].right
    File dense_vcf_shard = dense_vcf_info[i].left
    File dense_vcf_shard_idx = dense_vcf_info[i].right

    # Compute site-level metrics
    call QcTasks.CollectSiteMetrics {
      input:
        vcf = sites_vcf_shard,
        vcf_idx = sites_vcf_shard_idx,
        n_samples = ChooseTargetSamples.n_unrelated_samples,
        common_af_cutoff = common_af_cutoff,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Compute sample-level metrics
    call QcTasks.CollectSampleGenotypeMetrics {
      input:
        vcf = dense_vcf_shard,
        vcf_idx = dense_vcf_shard_idx,
        site_metrics = CollectSiteMetrics.all_sites,
        common_af_cutoff = common_af_cutoff,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Subset dense VCF to common variants 
    # (defined from all unrelated site freqs)
    call QcTasks.SubsetVcfByVids as CommonFilterDenseVcf {
      input:
        vcf = dense_vcf_shard,
        vcf_idx = dense_vcf_shard_idx,
        vids_list = CollectSiteMetrics.common_vids_list,
        bcftools_docker = bcftools_docker
    }
  }

  # Collect site-level metrics for dense subset if needed
  if ( has_twins || do_sample_bench ) {
    scatter ( i in range(length(dense_sites_vcf_info)) ) {
      call QcTasks.CollectSiteMetrics as DenseSiteMetrics {
        input:
          vcf = dense_sites_vcf_info[i].left,
          vcf_idx = dense_sites_vcf_info[i].right,
          n_samples = ChooseTargetSamples.n_target_samples,
          common_af_cutoff = common_af_cutoff,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Collapse all SNV sites
  Array[File] all_snv_beds = select_all(CollectSiteMetrics.snv_sites)
  if ( length(all_snv_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseAllSnvs {
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
    call QcTasks.ConcatTextFiles as CollapseCommonSnvs {
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
    call QcTasks.ConcatTextFiles as CollapseAllIndels {
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
    call QcTasks.ConcatTextFiles as CollapseCommonIndels {
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
    call QcTasks.ConcatTextFiles as CollapseAllSvs {
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
    call QcTasks.ConcatTextFiles as CollapseCommonSvs {
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

  # Collapse sample genotypes and split into one file per sample
  call QcTasks.ConcatGenotypeTsvs {
    input:
      tsvs = CollectSampleGenotypeMetrics.genotypes_tsv,
      output_prefix = output_prefix,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collapse all sites for dense subset
  Array[File] dense_snv_site_shards = select_all(select_first([DenseSiteMetrics.snv_sites, [empty_bed]]))
  Array[File] dense_indel_site_shards = select_all(select_first([DenseSiteMetrics.indel_sites, [empty_bed]]))
  Array[File] dense_sv_site_shards = select_all(select_first([DenseSiteMetrics.sv_sites, [empty_bed]]))
  Array[File] all_dense_sites_beds = flatten([dense_snv_site_shards, dense_indel_site_shards, dense_sv_site_shards])
  if ( length(all_dense_sites_beds) > 0 ) {
    call QcTasks.ConcatTextFiles as CollapseDenseSiteMetrics {
      input:
        shards = all_dense_sites_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = output_prefix + ".dense_all_sites.bed.gz",
        docker = bcftools_docker
    }
  }

  # Concatenate dense common VCFs
  # Necessary for LD comparisons (need all variants in a single file)
  call QcTasks.ConcatVcfs as ConcatCommonVcfs {
    input:
      vcfs = CommonFilterDenseVcf.subsetted_vcf,
      vcf_idxs = CommonFilterDenseVcf.subsetted_vcf_idx,
      out_prefix = output_prefix + ".dense.common",
      bcftools_concat_options = "--allow-overlaps",
      bcftools_docker = bcftools_docker
  }

  # Calculate peak LD between all pairs of variant classes
  call CalcLd as CalcCommonLd {
    input:
      vcf = ConcatCommonVcfs.merged_vcf,
      vcf_idx = ConcatCommonVcfs.merged_vcf_idx,
      common_snvs_bed = CollapseCommonSnvs.merged_file,
      common_indels_bed = CollapseCommonIndels.merged_file,
      common_svs_bed = CollapseCommonSvs.merged_file,
      out_prefix = output_prefix + ".dense.common",
      g2c_analysis_docker = g2c_analysis_docker
  }


  ################
  ### BENCHMARKING
  ################

  # Divide benchmarking shard quota evenly between site and sample benchmarking
  Int half_bench_shards = ceil(benchmarking_shards / 2)
  Int site_bench_shards = if do_sample_bench 
                          then ceil(half_bench_shards / (n_site_benchmark_datasets + 0.01))
                          else ceil(benchmarking_shards / (n_site_benchmark_datasets + 0.01))
  Int sample_bench_shards = if do_site_bench 
                            then ceil(half_bench_shards / (n_gt_benchmark_datasets + 0.01))
                            else ceil(benchmarking_shards / (n_gt_benchmark_datasets + 0.01))

  # Compute site-level benchmarking metrics
  if ( n_site_benchmark_datasets > 0 ) {
    scatter ( site_bench_idx in range(n_site_benchmark_datasets) ) {

      File sb_target_snv_bed = if length(snv_site_benchmark_beds) > site_bench_idx 
                               then select_all(snv_site_benchmark_beds)[site_bench_idx] 
                               else empty_bed
      File sb_target_indel_bed = if length(indel_site_benchmark_beds) > site_bench_idx
                                 then select_all(indel_site_benchmark_beds)[site_bench_idx]
                                 else empty_bed
      File sb_target_sv_bed = if length(sv_site_benchmark_beds) > site_bench_idx
                              then select_all(sv_site_benchmark_beds)[site_bench_idx]
                              else empty_bed
      String sb_target_prefix = if defined(site_benchmark_dataset_names) 
                                then select_all(site_benchmark_dataset_names)[site_bench_idx] 
                                else ""

      call BenchSites.BenchmarkSites {
        input:
          source_snv_bed = CollapseAllSnvs.merged_file,
          source_indel_bed = CollapseAllIndels.merged_file,
          source_sv_bed = CollapseAllSvs.merged_file,
          source_prefix = output_prefix,
          target_snv_bed = sb_target_snv_bed,
          target_indel_bed = sb_target_indel_bed,
          target_sv_bed = sb_target_sv_bed,
          target_prefix = sb_target_prefix,
          eval_interval_beds = select_first([benchmark_interval_beds]),
          eval_interval_bed_names = select_first([benchmark_interval_bed_names]),
          genome_file = select_first([genome_file]),
          total_shards = site_bench_shards,
          common_af_cutoff = common_af_cutoff,
          bcftools_docker = bcftools_docker,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Compute sample-level benchmarking metrics
  if ( n_gt_benchmark_datasets > 0 ) {
    scatter ( gt_bench_idx in range(n_gt_benchmark_datasets) ) {

      Array[File] gb_target_vcfs = if length(sample_benchmark_id_maps) > gt_bench_idx
                                   then select_all(sample_benchmark_vcfs[gt_bench_idx])
                                   else [empty_vcf]
      Array[File] gb_target_vcf_idxs = if length(sample_benchmark_id_maps) > gt_bench_idx
                                       then select_all(sample_benchmark_vcf_idxs[gt_bench_idx])
                                       else [empty_vcf_idx]
      Array[File] gb_target_id_maps = if length(sample_benchmark_id_maps) > gt_bench_idx
                                      then select_all(sample_benchmark_id_maps[gt_bench_idx])
                                      else [MakeEmptyFile.empty_file]
      String gb_target_prefix = if length(sample_benchmark_dataset_names) > gt_bench_idx
                                then select_all(sample_benchmark_dataset_names)[gt_bench_idx]
                                else ""

      call BenchSamples.BenchmarkSamples {
        input:
          source_all_sites_bed = CollapseDenseSiteMetrics.merged_file,
          source_snv_beds = select_first([DenseSiteMetrics.snv_sites, empty_bed]),
          source_indel_beds = select_first([DenseSiteMetrics.indel_sites, empty_bed]),
          source_sv_beds = select_first([DenseSiteMetrics.sv_sites, empty_bed]),
          source_gt_tarball = ConcatGenotypeTsvs.genotypes_tarball,
          source_prefix = output_prefix,
          target_vcfs = gb_target_vcfs,
          target_vcf_idxs = gb_target_vcf_idxs,
          target_prefix = gb_target_prefix,
          id_map_tsvs = gb_target_id_maps,
          eval_interval_beds = select_first([benchmark_interval_beds]),
          eval_interval_bed_names = select_first([benchmark_interval_bed_names]),
          genome_file = select_first([genome_file]),
          total_shards = sample_bench_shards,
          min_samples_per_shard = min_samples_per_bench_shard,
          common_af_cutoff = common_af_cutoff,
          bcftools_docker = bcftools_docker,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Perform genotype benchmarking between twins/technical replicates, if available
  if ( has_twins ) {
    call BenchTwins.BenchmarkTwins {
      input:
        all_sites_bed = select_first([CollapseDenseSiteMetrics.merged_file, empty_bed]),
        gt_tarball = ConcatGenotypeTsvs.genotypes_tarball,
        twins_tsv = select_first([CleanTwins.complete_twins_tsv, MakeEmptyFile.empty_file]),
        output_prefix = "~{output_prefix}.twins_techreps",
        eval_interval_beds = select_first([benchmark_interval_beds]),
        eval_interval_bed_names = select_first([benchmark_interval_bed_names]),
        genome_file = select_first([genome_file]),
        common_af_cutoff = common_af_cutoff,
        bcftools_docker = bcftools_docker,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # Perform trio-based Mendelian benchmarking, if available
  if ( has_trios ) {
    call BenchTrios.BenchmarkTrios {
      input:
        vcfs = dense_vcf_shards,
        vcf_idxs = dense_vcf_shard_idxs,
        sites_bed = select_first([CollapseDenseSiteMetrics.merged_file]),
        trios_fam = select_first([CleanFam.trios_fam]),
        trios_samples_list = select_first([CleanFam.trio_samples_list]),
        eval_interval_beds = select_first([benchmark_interval_beds]),
        eval_interval_bed_names = select_first([benchmark_interval_bed_names]),
        genome_file = select_first([genome_file]),
        output_prefix = "~{output_prefix}.trios",
        common_af_cutoff = common_af_cutoff,
        bcftools_docker = bcftools_docker,
        g2c_analysis_docker = g2c_analysis_docker
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

  # Collapse genotype distributions
  call QcTasks.SumCompressedDistribs as SumGenotypeDistribs {
    input:
      distrib_tsvs = select_all(CollectSampleGenotypeMetrics.compressed_gt_distrib),
      n_key_columns = 4,
      out_prefix = output_prefix + ".genotype_distribution",
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
    File genotype_distrib = SumGenotypeDistribs.merged_distrib
    File ld_stats = CalcCommonLd.ld_stats

    Array[Array[File]]? site_benchmark_common_snv_ppv_beds = BenchmarkSites.common_snv_ppv_beds
    Array[Array[File]]? site_benchmark_common_snv_sens_beds = BenchmarkSites.common_snv_sens_beds
    Array[Array[File]]? site_benchmark_common_indel_ppv_beds = BenchmarkSites.common_indel_ppv_beds
    Array[Array[File]]? site_benchmark_common_indel_sens_beds = BenchmarkSites.common_indel_sens_beds
    Array[Array[File]]? site_benchmark_common_sv_ppv_beds = BenchmarkSites.common_sv_ppv_beds
    Array[Array[File]]? site_benchmark_common_sv_sens_beds = BenchmarkSites.common_sv_sens_beds
    Array[Array[File]]? site_benchmark_ppv_by_sizes = BenchmarkSites.ppv_by_sizes
    Array[Array[File]]? site_benchmark_sensitivity_by_sizes = BenchmarkSites.sensitivity_by_sizes
    Array[Array[File]]? site_benchmark_ppv_by_freqs = BenchmarkSites.ppv_by_freqs
    Array[Array[File]]? site_benchmark_sensitivity_by_freqs = BenchmarkSites.sensitivity_by_freqs
    Array[Array[File]]? site_benchmark_false_positive_runs = BenchmarkSites.fp_runs
    Array[Array[File]]? site_benchmark_false_negative_runs = BenchmarkSites.fn_runs

    Array[Array[File]]? sample_benchmark_ppv_distribs = BenchmarkSamples.ppv_distribs
    Array[Array[File]]? sample_benchmark_sensitivity_distribs = BenchmarkSamples.sensitivity_distribs
    Array[File]? twin_genotype_benchmark_distribs = BenchmarkTwins.gt_benchmarking_distribs
    Array[File]? trio_mendelian_violation_distribs = BenchmarkTrios.trio_benchmarking_distribs
  }
}


# Calculate peak LD R2 statistics per variant between all pairs of variants
task CalcLd {
  input {
    File vcf
    File vcf_idx
    
    File? common_snvs_bed
    File? common_indels_bed
    File? common_svs_bed

    Int plink_ld_window_kb = 1000
    Float plink_ld_window_min_r2 = 0.05
    
    String out_prefix
    String g2c_analysis_docker
  }

  Boolean has_snvs = defined(common_snvs_bed)
  Boolean has_indels = defined(common_indels_bed)
  Boolean has_svs = defined(common_svs_bed)

  Int disk_gb_auto = ceil(20 * size(vcf, "GB")) + 10
  Int disk_gb = if disk_gb_auto > 1000 then 1000 else disk_gb_auto

  command <<<
    set -eu -o pipefail

    # Write list of variant IDs present in VCF
    bcftools query -f '%ID\n' ~{vcf} > elig_vids.list

    # Define lists of common variants for each class
    if ~{has_snvs}; then
      zcat ~{default="" common_snvs_bed} \
      | cut -f4 | fgrep -xf elig_vids.list \
      > snvs.list || true
    fi
    if ~{has_indels}; then
      zcat ~{default="" common_indels_bed} \
      | cut -f4 | fgrep -xf elig_vids.list \
      > indels.list || true
    fi
    if ~{has_svs}; then
      zcat ~{default="" common_svs_bed} \
      | cut -f4 | fgrep -xf elig_vids.list \
      > svs.list || true
    fi

    # Compute pairwise LD for all pairs of variants
    plink2 \
      --r2-unphased 'yes-really' \
      --ld-window-kb ~{plink_ld_window_kb} \
      --ld-window-r2 ~{plink_ld_window_min_r2} \
      --vcf ~{vcf} \
      --out ~{out_prefix}

    # Simulate synthetic LD "noise" for each variant
    # Note that this R script appends to --out-tsv, so all 
    # variant classes can be written to the same file without issue
    n_samples=$( bcftools query -l ~{vcf} | wc -l )
    if ~{has_snvs}; then
      /opt/pancan_germline_wgs/scripts/qc/vcf_qc/null_ld_sim.R \
        --bed ~{default="" common_snvs_bed} \
        --sample-size $n_samples \
        --max-r2 ~{plink_ld_window_min_r2} \
        --out-tsv null_r2.tsv
    fi
    if ~{has_indels}; then
      /opt/pancan_germline_wgs/scripts/qc/vcf_qc/null_ld_sim.R \
        --bed ~{default="" common_indels_bed} \
        --sample-size $n_samples \
        --max-r2 ~{plink_ld_window_min_r2} \
        --out-tsv null_r2.tsv
    fi
    if ~{has_svs}; then
      /opt/pancan_germline_wgs/scripts/qc/vcf_qc/null_ld_sim.R \
        --bed ~{default="" common_svs_bed} \
        --sample-size $n_samples \
        --max-r2 ~{plink_ld_window_min_r2} \
        --out-tsv null_r2.tsv
    fi

    # Iterate over pairs of variant classes and 
    # report the strongest R2 per variant
    for vc1 in snv indel sv; do
      if [ -s $vc1.list ]; then
        for vc2 in snv indel sv; do
          if [ -s $vc2.list ]; then
            fgrep -wf $vc1.list \
              ~{out_prefix}.ld.vcor \
            | fgrep -wf $vc2.list \
            | awk -v FS="\t" -v OFS="\t" \
              '{ if ($3!=$6) print $3, $7"\n"$6, $7 }' \
            | cat - null_r2.tsv \
            | fgrep -wf $vc1.list \
            | sort -nrk2,2 -k1,1V \
            | awk '!seen[$1]++' \
            | awk -v OFS="\t" -v vc2=$vc2 '{ print $1, vc2, $2 }'
          fi
        done
      fi
    done \
    | sort -Vk1,1 -k2,2V \
    | cat <( echo -e "#vid\tother_vc\tld_r2" ) - | gzip -c \
    > ~{out_prefix}.peak_ld_by_vc.tsv.gz
  >>>

  output {
    File ld_stats = "~{out_prefix}.peak_ld_by_vc.tsv.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }  
}


# Define lists of "target" samples to use for site-level and sample-based analyses
task ChooseTargetSamples {
  input {
    File all_samples_list
    File? probands_list
    File? duplicate_samples_list
    String out_prefix
    File? sample_priority_tsv
    Int n_samples
    String g2c_analysis_docker
  }

  String site_excl_outfile = out_prefix + ".exclude_for_site_analysis.samples.list"
  String target_outfile = out_prefix + ".target.samples.list"

  command <<<
    set -eu -o pipefail

    # Debugging
    echo -e "Contents of working directory:"
    find ./

    # If probands or duplicates are provided, exclude them from sample universe
    echo -e "THIS_SAMPLE_SHOULD_NEVER_HIT" > site.exclude.samples.list
    if ~{defined(probands_list)}; then
      cat ~{default="" probands_list} >> site.exclude.samples.list || true
    fi
    if ~{defined(duplicate_samples_list)}; then
      cat ~{default="" duplicate_samples_list} >> site.exclude.samples.list || true
    fi
    sort -V site.exclude.samples.list | uniq > ~{site_excl_outfile}

    # Produce list of unrelated samples (for counting effective sample size only)
    fgrep \
      -xvf site.exclude.samples.list \
      ~{all_samples_list} \
    | sort -V | uniq \
    > unrelated.samples.list

    # Select samples to include for genotype-level analyses
    cmd="/opt/pancan_germline_wgs/scripts/qc/vcf_qc/select_qc_samples.R"
    cmd="$cmd --all-samples-list \"~{all_samples_list}\" -N ~{n_samples}"
    cmd="$cmd --out-list \"~{target_outfile}\""
    if ~{defined(sample_priority_tsv)}; then
      cmd="$cmd --priority-tsv \"~{default="" sample_priority_tsv}\""
    fi
    echo -e "Now selecting samples with this command:\n$cmd"
    eval $cmd
  >>>

  output {
    File target_samples = target_outfile
    File site_exclude_samples = site_excl_outfile
    Int n_unrelated_samples = length(read_lines("unrelated.samples.list"))
    Int n_target_samples = length(read_lines(target_outfile))
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk 20 HDD"
    preemptible: 3
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
  String trio_samples_out_fname = basename(fam_file, ".fam") + ".trio_samples.list"
  
  command <<<
    set -eu -o pipefail

    # Clean .fam file to duos or trios present in all_samples_list
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/subset_fam.R \
      --in-fam ~{fam_file} \
      --all-samples ~{all_samples_list} \
      --out-fam duos_and_trios.fam

    # Write list of all probands (can exclude these samples for better AF estimates)
    fgrep -v "#" duos_and_trios.fam | cut -f2 > "~{proband_out_fname}"

    # Further restrict .fam to complete trios
    awk -v FS="\t" -v OFS="\t" \
      '{ if ($2!=0 && $3!=0 && $4!=0) print }' \
      duos_and_trios.fam \
    > "~{fam_out_fname}"

    # Write list of all samples in complete trios
    awk -v OFS="\n" '{ print $2, $3, $4 }' "~{fam_out_fname}" \
    | sort -V | uniq \
    > "~{trio_samples_out_fname}"
  >>>

  output {
    File trios_fam = "~{fam_out_fname}"
    Int n_trios = length(read_lines("~{fam_out_fname}"))
    File probands_list = "~{proband_out_fname}"
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


# Clean input .tsv file to restrict to complete twin pairs present in input VCF
task CleanTwins {
  input {
    File twins_tsv
    File all_samples_list
    String docker
  }

  String twins_out_fname = basename(twins_tsv, ".tsv") + ".cleaned_twins.tsv"
  String duplicate_out_fname = basename(twins_tsv, ".tsv") + ".duplicate_ids.list"
  String twin_samples_out_fname = basename(twins_tsv, ".tsv") + ".twin_pair_samples.list"
  
  command <<<
    set -eu -o pipefail

    # Coerce twins .tsv to fake .fam file for cleaning
    awk -v FS="\t" -v OFS="\t" \
      '{ print NR, $1, $2, "0", "0", "0" }' \
      ~{twins_tsv} \
    > twins.fam

    # Clean .fam file to duos or trios present in all_samples_list
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/subset_fam.R \
      --in-fam twins.fam \
      --all-samples ~{all_samples_list} \
      --out-fam cleaned_twins.fam

    # Write list of complete twin pairs
    awk -v FS="\t" -v OFS="\t" \
      '{ if ($2!=0 && $3!=0) print $2, $3 }' \
      cleaned_twins.fam \
    > "~{twins_out_fname}"

    # Write list of all samples involved in complete trios
    awk -v OFS="\n" '{ print $1, $2 }' "~{twins_out_fname}" \
    | sort -V | uniq \
    > "~{twin_samples_out_fname}"

    # Write list of one member from each pair (can exclude these samples for better AF estimates)
    cut -f1 "~{twins_out_fname}" | sort | uniq > "~{duplicate_out_fname}"
  >>>

  output {
    File complete_twins_tsv = "~{twins_out_fname}"
    Int n_twins = length(read_lines("~{twins_out_fname}"))
    File duplicate_samples_list = "~{duplicate_out_fname}"
    File twin_samples_list = "~{twin_samples_out_fname}"
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
    File? site_exclude_samples
    Boolean has_mcnvs = false
    String extra_commands = ""
    File supp_vcf_header
    
    String out_prefix
    
    Int? disk_gb
    Float mem_gb = 4
    Int n_cpu = 2
    
    String g2c_analysis_docker
  }

  String no_rel_cmd = if defined(site_exclude_samples) then "--force-samples --samples-file ^" + basename(select_first([site_exclude_samples])) else ""
  String sites_outfile = out_prefix + ".sites.vcf.gz"
  String dense_outfile = out_prefix + ".dense_subset.vcf.gz"
  String dense_sites_outfile = out_prefix + ".dense_subset.sites.vcf.gz"

  String mcnv_anno = if has_mcnvs then "| /opt/pancan_germline_wgs/scripts/gatksv_helpers/annotate_mcnv_freqs.py - -" else ""

  Int default_disk_gb = ceil(4 * size(vcf, "GB")) + 10
  Int hdd_gb = select_first([disk_gb, default_disk_gb])

  command <<<
    set -eu -o pipefail

    # Move site exclude samples to working directory, if optioned
    if ~{defined(site_exclude_samples)}; then
      mv ~{select_first([site_exclude_samples])} ./
    fi

    # Generate sites-only VCF of unrelated samples
    bcftools view ~{no_rel_cmd} ~{vcf} \
    | bcftools annotate -h ~{supp_vcf_header} --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
    | bcftools +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Het,AC_Hom,HWE \
    ~{mcnv_anno} \
    ~{extra_commands} \
    | bcftools annotate -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE,^FILTER/PASS,FILTER/MULTIALLELIC" \
    | bcftools view -G --threads 2 \
      --include 'INFO/AC > 0 | FILTER="MULTIALLELIC"' \
      -Oz -o ~{sites_outfile}
    tabix -p vcf -f ~{sites_outfile}

    # Generate dense VCF of only target samples
    # Note that this command also reassigns all variant IDs to G2C QC standards
    # This is required for compatability with downstream LD calculations
    bcftools view --samples-file ~{target_samples} --force-samples ~{vcf} \
    | bcftools annotate -h ~{supp_vcf_header} --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
    | bcftools view --include 'INFO/AC > 0 | FILTER="MULTIALLELIC"' \
    | bcftools +fill-tags -- -t AN,AC,AF,AC_Hemi,AC_Het,AC_Hom,HWE \
    ~{mcnv_anno} \
    ~{extra_commands} \
    | bcftools annotate -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE,^FILTER/PASS,FILTER/MULTIALLELIC,^FORMAT/GT,FORMAT/RD_CN" \
    | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/set_g2c_qc_variant_ids.py \
      --vcf-out ~{dense_outfile}
    tabix -p vcf -f ~{dense_outfile}

    # Make a sites-only equivalent of the dense VCF for downstream benchmarking
    bcftools view -G ~{dense_outfile} -Oz -o ~{dense_sites_outfile}
    tabix -p vcf -f ~{dense_sites_outfile}

    # For convenience, also make an empty VCF and corresponding index with matching header
    tabix --only-header ~{dense_outfile} | bgzip -c > empty.vcf.gz
    tabix -p vcf -f empty.vcf.gz
  >>>

  output {
    File sites_vcf = sites_outfile
    File sites_vcf_idx = "~{sites_outfile}.tbi"
    File dense_vcf = dense_outfile
    File dense_vcf_idx = "~{dense_outfile}.tbi"
    File dense_sites_vcf = dense_sites_outfile
    File dense_sites_vcf_idx = "~{dense_sites_outfile}.tbi"
    File empty_vcf = "empty.vcf.gz"
    File empty_vcf_idx = "empty.vcf.gz.tbi"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk ~{hdd_gb} HDD"
    preemptible: 3
  }
}

