# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# WDL to run first half of GATK "biggest practices" / "Gnarly" joint genotyping of SNVs/indels from gVCF input

# Unlike the parent WARP workflow, this task is designed to be executed on smaller subsets of sharded intervals
# E.g., per chromosome (or even smaller regions)
# These individual shards can then be merged prior to VQSR, fingerprinting, etc
# See GnarlyJointGenotypingPart2.wdl for those subsequent steps

# Based on GATK WARP pipeline:
# https://github.com/broadinstitute/warp/blob/develop/pipelines/broad/dna_seq/germline/joint_genotyping/JointGenotyping.wdl

# See also:
# https://gatk.broadinstitute.org/hc/en-us/articles/16957867036315-Introducing-GATK-Biggest-Practices-for-Joint-Calling-Supersized-Cohorts


version 1.0


import "https://raw.githubusercontent.com/broadinstitute/warp/develop/tasks/broad/JointGenotypingTasks.wdl" as Tasks
import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/refs/heads/gatkhc/wdl/Utilities.wdl" as G2CUtils


workflow GnarlyJointGenotypingPart1 {
  input {
    File unpadded_intervals_file

    String callset_name
    File sample_name_map
    Int import_gvcfs_batch_size = 50

    File ref_fasta
    File ref_fasta_index
    File ref_dict

    File dbsnp_vcf

    Boolean make_hard_filtered_sites = true

    Int small_disk = 100
    Int import_gvcfs_disk_gb = 200
    Int make_sites_disk_gb = 200
    Int large_disk = 1000

    Int? top_level_scatter_count
    Boolean intervals_already_split = false
    Float unbounded_scatter_count_scale_factor = 0.15
    Int gnarly_scatter_count = 10
  }

  Array[Array[String]] sample_name_map_lines = read_tsv(sample_name_map)
  Int num_gvcfs = length(sample_name_map_lines)

  Int unbounded_scatter_count = select_first([top_level_scatter_count, round(unbounded_scatter_count_scale_factor * num_gvcfs)])
  Int scatter_count = if unbounded_scatter_count > 2 then unbounded_scatter_count else 2


  call Tasks.CheckSamplesUnique {
    input:
      sample_name_map = sample_name_map
  }

  # We deviate from the WARP protocol by enabling the user to provide custom pre-split intervals
  # This requires a custom task to strictly split the interval list into single-interval shards
  if ( intervals_already_split ) {
    call G2CUtils.SplitIntervalList as G2CSplitIntervals {
      input:
        interval_list = unpadded_intervals_file
    }
  }
  if ( !intervals_already_split ) {
    call Tasks.SplitIntervalList as GatkSplitIntervals {
      input:
        interval_list = unpadded_intervals_file,
        scatter_count = scatter_count,
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        disk_size_gb = small_disk,
        sample_names_unique_done = CheckSamplesUnique.samples_unique
    }
  }
  
  Array[File] unpadded_intervals = select_first([G2CSplitIntervals.output_intervals, GatkSplitIntervals.output_intervals])

  scatter (idx in range(length(unpadded_intervals))) {
    call Tasks.ImportGVCFs {
      input:
        sample_name_map = sample_name_map,
        interval = unpadded_intervals[idx],
        ref_fasta = ref_fasta,
        ref_fasta_index = ref_fasta_index,
        ref_dict = ref_dict,
        workspace_dir_name = "genomicsdb",
        disk_size_gb = import_gvcfs_disk_gb,
        batch_size = import_gvcfs_batch_size
    }

    if ( gnarly_scatter_count > 1 ) {
      call Tasks.SplitIntervalList as GnarlyIntervalScatterDude {
        input:
          interval_list = unpadded_intervals[idx],
          scatter_count = gnarly_scatter_count,
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          disk_size_gb = small_disk,
          sample_names_unique_done = CheckSamplesUnique.samples_unique
      }
    }

    Array[File] gnarly_intervals = select_first([GnarlyIntervalScatterDude.output_intervals, [unpadded_intervals[idx]]])

    scatter (gnarly_idx in range(length(gnarly_intervals))) {
      call Tasks.GnarlyGenotyper {
        input:
          workspace_tar = ImportGVCFs.output_genomicsdb,
          interval = gnarly_intervals[gnarly_idx],
          output_vcf_filename = callset_name + "." + idx + "." + gnarly_idx + ".vcf.gz",
          ref_fasta = ref_fasta,
          ref_fasta_index = ref_fasta_index,
          ref_dict = ref_dict,
          dbsnp_vcf = dbsnp_vcf
      }
    }

    Array[File] gnarly_vcfs = GnarlyGenotyper.output_vcf
    Array[File] gnarly_vcf_idxs = GnarlyGenotyper.output_vcf_index

    if ( length(gnarly_vcfs) > 1 ) {
      call Tasks.GatherVcfs as TotallyRadicalGatherVcfs {
        input:
          input_vcfs = gnarly_vcfs,
          output_vcf_name = callset_name + "." + idx + ".gnarly.vcf.gz",
          disk_size_gb = large_disk
      }
    }

    File genotyped_vcf = select_first([TotallyRadicalGatherVcfs.output_vcf, gnarly_vcfs[0]])
    File genotyped_vcf_index = select_first([TotallyRadicalGatherVcfs.output_vcf_index, gnarly_vcf_idxs[0]])

    if ( make_hard_filtered_sites ){
      call Tasks.HardFilterAndMakeSitesOnlyVcf {
        input:
          vcf = genotyped_vcf,
          vcf_index = genotyped_vcf_index,
          excess_het_threshold = 54.69,
          variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
          sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
          disk_size_gb = make_sites_disk_gb
      }
    }

    File jg_vcf_shard = select_first([HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf, genotyped_vcf])
    File jg_vcf_shard_idx = select_first([HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index, genotyped_vcf_index])
  }

  output {
    Array[File] joint_genotyped_vcf = jg_vcf_shard
    Array[File] joint_genotyped_vcf_index = jg_vcf_shard_idx
    Array[File?] sites_only_vcfs = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf
    Array[File?] sites_only_vcfs_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index
  }
}

