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


workflow GnarlyJointGenotypingPart2 {
  input {
    # Below are the inputs that deviate from the WARP joint genotyping pipeline
    # Note that these are provided as nested arrays because step 1 is intended to 
    # be parallelized (e.g., per chromosome). These are flattened below.
    Array[Array[File]] sites_only_vcfs_nested
    Array[Array[File]] sites_only_vcfs_index_nested
    Array[Array[File]] variant_filtered_vcfs_nested
    Array[Array[File]] variant_filtered_vcfs_index_nested

    # All other inputs below this point are identical to the public WARP pipeline
    String callset_name
    File sample_name_map

    File ref_dict

    Int small_disk = 100
    Int medium_disk = 200
    Int large_disk = 1000
    Int huge_disk = 2000

    Array[String] snp_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.8", "99.6", "99.5", "99.4", "99.3", "99.0", "98.0", "97.0", "90.0" ]
    Array[String] snp_recalibration_annotation_values = ["QD", "MQRankSum", "ReadPosRankSum", "FS", "MQ", "SOR", "DP"]
    Array[String] indel_recalibration_tranche_values = ["100.0", "99.95", "99.9", "99.5", "99.0", "97.0", "96.0", "95.0", "94.0", "93.5", "93.0", "92.0", "91.0", "90.0"]
    Array[String] indel_recalibration_annotation_values = ["FS", "ReadPosRankSum", "MQRankSum", "QD", "SOR", "DP"]

    File eval_interval_list
    File hapmap_resource_vcf
    File hapmap_resource_vcf_index
    File omni_resource_vcf
    File omni_resource_vcf_index
    File one_thousand_genomes_resource_vcf
    File one_thousand_genomes_resource_vcf_index
    File mills_resource_vcf
    File mills_resource_vcf_index
    File axiomPoly_resource_vcf
    File axiomPoly_resource_vcf_index
    File dbsnp_resource_vcf
    File dbsnp_resource_vcf_index

    Float snp_filter_level = 99.7
    Float indel_filter_level = 99
    Int snp_vqsr_downsampleFactor = 10
  }

  Int num_gvcfs = length(read_tsv(sample_name_map))

  Boolean is_small_callset = num_gvcfs <= 1000

  Array[File] sites_only_vcfs = flatten(sites_only_vcfs_nested)
  Array[File] sites_only_vcfs_index = flatten(sites_only_vcfs_index_nested)
  Array[File] variant_filtered_vcfs = flatten(variant_filtered_vcfs_nested)
  Array[File] variant_filtered_vcfs_index = flatten(variant_filtered_vcfs_index_nested)

  call Tasks.GatherVcfs as SitesOnlyGatherVcf {
    input:
      input_vcfs = sites_only_vcfs,
      output_vcf_name = callset_name + ".sites_only.vcf.gz",
      disk_size_gb = medium_disk
  }

  call Tasks.IndelsVariantRecalibrator {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".indels.recal",
      tranches_filename = callset_name + ".indels.tranches",
      recalibration_tranche_values = indel_recalibration_tranche_values,
      recalibration_annotation_values = indel_recalibration_annotation_values,
      mills_resource_vcf = mills_resource_vcf,
      mills_resource_vcf_index = mills_resource_vcf_index,
      axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      use_allele_specific_annotations = false,
      disk_size_gb = small_disk
  }

  call Tasks.SNPsVariantRecalibratorCreateModel {
    input:
      sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      recalibration_filename = callset_name + ".snps.recal",
      tranches_filename = callset_name + ".snps.tranches",
      recalibration_tranche_values = snp_recalibration_tranche_values,
      recalibration_annotation_values = snp_recalibration_annotation_values,
      downsampleFactor = snp_vqsr_downsampleFactor,
      model_report_filename = callset_name + ".snps.model.report",
      hapmap_resource_vcf = hapmap_resource_vcf,
      hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      omni_resource_vcf = omni_resource_vcf,
      omni_resource_vcf_index = omni_resource_vcf_index,
      one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      dbsnp_resource_vcf = dbsnp_resource_vcf,
      dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
      use_allele_specific_annotations = false,
      disk_size_gb = small_disk
  }

  scatter (idx in range(length(sites_only_vcfs))) {
    call Tasks.SNPsVariantRecalibrator as SNPsVariantRecalibratorScattered {
      input:
        sites_only_variant_filtered_vcf = sites_only_vcfs[idx],
        sites_only_variant_filtered_vcf_index = sites_only_vcfs_index[idx],
        recalibration_filename = callset_name + ".snps." + idx + ".recal",
        tranches_filename = callset_name + ".snps." + idx + ".tranches",
        recalibration_tranche_values = snp_recalibration_tranche_values,
        recalibration_annotation_values = snp_recalibration_annotation_values,
        model_report = SNPsVariantRecalibratorCreateModel.model_report,
        hapmap_resource_vcf = hapmap_resource_vcf,
        hapmap_resource_vcf_index = hapmap_resource_vcf_index,
        omni_resource_vcf = omni_resource_vcf,
        omni_resource_vcf_index = omni_resource_vcf_index,
        one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
        one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
        dbsnp_resource_vcf = dbsnp_resource_vcf,
        dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
        use_allele_specific_annotations = false,
        disk_size_gb = small_disk
      }
    }

  call Tasks.GatherTranches as SNPGatherTranches {
    input:
      tranches = SNPsVariantRecalibratorScattered.tranches,
      output_filename = callset_name + ".snps.gathered.tranches",
      mode = "SNP",
      disk_size_gb = small_disk
  }

  scatter (idx in range(length(variant_filtered_vcfs))) {
    call Tasks.ApplyRecalibration {
      input:
        recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
        input_vcf = variant_filtered_vcfs[idx],
        input_vcf_index = variant_filtered_vcfs_index[idx],
        indels_recalibration = IndelsVariantRecalibrator.recalibration,
        indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
        indels_tranches = IndelsVariantRecalibrator.tranches,
        snps_recalibration = SNPsVariantRecalibratorScattered.recalibration[idx],
        snps_recalibration_index = SNPsVariantRecalibratorScattered.recalibration_index[idx],
        snps_tranches = SNPGatherTranches.tranches_file,
        indel_filter_level = indel_filter_level,
        snp_filter_level = snp_filter_level,
        use_allele_specific_annotations = false,
        disk_size_gb = medium_disk
    }

    # For large callsets we need to collect metrics from the shards and gather them later.
    if (!is_small_callset) {
      call Tasks.CollectVariantCallingMetrics as CollectMetricsSharded {
        input:
          input_vcf = ApplyRecalibration.recalibrated_vcf,
          input_vcf_index = ApplyRecalibration.recalibrated_vcf_index,
          metrics_filename_prefix = callset_name + "." + idx,
          dbsnp_vcf = dbsnp_resource_vcf,
          dbsnp_vcf_index = dbsnp_resource_vcf_index,
          interval_list = eval_interval_list,
          ref_dict = ref_dict,
          disk_size_gb = medium_disk
      }
    }
  }

  # For small callsets we can gather the VCF shards and then collect metrics on it.
  if (is_small_callset) {
    call Tasks.GatherVcfs as FinalGatherVcf {
      input:
        input_vcfs = ApplyRecalibration.recalibrated_vcf,
        output_vcf_name = callset_name + ".vcf.gz",
        disk_size_gb = huge_disk
    }

    call Tasks.CollectVariantCallingMetrics as CollectMetricsOnFullVcf {
      input:
        input_vcf = FinalGatherVcf.output_vcf,
        input_vcf_index = FinalGatherVcf.output_vcf_index,
        metrics_filename_prefix = callset_name,
        dbsnp_vcf = dbsnp_resource_vcf,
        dbsnp_vcf_index = dbsnp_resource_vcf_index,
        interval_list = eval_interval_list,
        ref_dict = ref_dict,
        disk_size_gb = large_disk
    }
  }

  if (!is_small_callset) {
    call Tasks.GatherVariantCallingMetrics {
      input:
        input_details = select_all(CollectMetricsSharded.detail_metrics_file),
        input_summaries = select_all(CollectMetricsSharded.summary_metrics_file),
        output_prefix = callset_name,
        disk_size_gb = medium_disk
    }
  }

  # Get the metrics from either code path
  File output_detail_metrics_file = select_first([CollectMetricsOnFullVcf.detail_metrics_file, GatherVariantCallingMetrics.detail_metrics_file])
  File output_summary_metrics_file = select_first([CollectMetricsOnFullVcf.summary_metrics_file, GatherVariantCallingMetrics.summary_metrics_file])

  # Get the VCFs from either code path
  Array[File?] output_vcf_files = if defined(FinalGatherVcf.output_vcf) then [FinalGatherVcf.output_vcf] else ApplyRecalibration.recalibrated_vcf
  Array[File?] output_vcf_index_files = if defined(FinalGatherVcf.output_vcf_index) then [FinalGatherVcf.output_vcf_index] else ApplyRecalibration.recalibrated_vcf_index

  output {
    # Metrics from either the small or large callset
    File detail_metrics_file = output_detail_metrics_file
    File summary_metrics_file = output_summary_metrics_file

    # Outputs from the small callset path through the wdl.
    Array[File] output_vcfs = select_all(output_vcf_files)
    Array[File] output_vcf_indices = select_all(output_vcf_index_files)
  }
  meta {
    allowNestedInputs: true
  }
}
