# Drop-in replacement for GATK-SV v1.0 module 08
# Simplified for DFCI G2C

version 1.0

import "FilterOutlierSamples.wdl" as filter_outliers
import "Structs.wdl"
import "Utils.wdl" as util
import "TasksMakeCohortVcf.wdl" as tasks_mcv

# Workflow to filter a pre-defined list of outliers from VCFs as part of FilterBatch after FilterBatchSites
workflow FilterBatchSamples {
  input {
    String batch
    File? manta_vcf
    File? wham_vcf
    File? melt_vcf
    File? depth_vcf
    File outlier_sample_list
    String sv_base_mini_docker
    String linux_docker
    RuntimeAttr? runtime_attr_subset_vcf
    RuntimeAttr? runtime_attr_ids_from_vcf
    RuntimeAttr? runtime_attr_filter_samples
    RuntimeAttr? runtime_attr_count_svs
    RuntimeAttr? runtime_attr_merge_pesr_vcfs
  }

  Array[File?] vcfs = [manta_vcf, wham_vcf, melt_vcf, depth_vcf]
  Array[String] algorithms = ["manta", "wham", "melt", "depth"]  # fixed algorithms to enable File outputs to be determined
  Int num_algorithms = length(algorithms)

  scatter (i in range(num_algorithms)) {
    if (defined(vcfs[i])) {
      call util.SubsetVcfBySamplesList {
        input:
          vcf = select_first([vcfs[i]]),
          list_of_samples = outlier_sample_list,
          outfile_name = "${batch}.${algorithms[i]}.outliers_removed.vcf.gz",
          remove_samples = true,
          sv_base_mini_docker = sv_base_mini_docker,
          runtime_attr_override = runtime_attr_subset_vcf
      }
    }
  }
  
  call util.GetSampleIdsFromVcf {
    input:
      vcf = select_first(vcfs),
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_ids_from_vcf
  }

  # Write new list of samples without outliers
  call filter_outliers.FilterSampleList {
    input:
      original_samples = GetSampleIdsFromVcf.out_array,
      outlier_samples = read_lines(outlier_sample_list),
      batch = batch,
      linux_docker = linux_docker,
      runtime_attr_override = runtime_attr_filter_samples
  }

  Array[File] pesr_vcfs_no_outliers = select_all([SubsetVcfBySamplesList.vcf_subset[0], SubsetVcfBySamplesList.vcf_subset[1], SubsetVcfBySamplesList.vcf_subset[2]])
  Array[File] pesr_vcfs_no_outliers_index = select_all([SubsetVcfBySamplesList.vcf_subset_index[0], SubsetVcfBySamplesList.vcf_subset_index[1], SubsetVcfBySamplesList.vcf_subset_index[2]])
  call tasks_mcv.ConcatVcfs as MergePesrVcfs {
    input:
      vcfs = pesr_vcfs_no_outliers,
      vcfs_idx = pesr_vcfs_no_outliers_index,
      allow_overlaps = true,
      outfile_prefix = "~{batch}.filtered_pesr_merged",
      sv_base_mini_docker = sv_base_mini_docker,
      runtime_attr_override = runtime_attr_merge_pesr_vcfs
  }

  output {
    File? outlier_filtered_manta_vcf = SubsetVcfBySamplesList.vcf_subset[0]
    File? outlier_filtered_wham_vcf = SubsetVcfBySamplesList.vcf_subset[1]
    File? outlier_filtered_melt_vcf = SubsetVcfBySamplesList.vcf_subset[2]
    File? outlier_filtered_depth_vcf = SubsetVcfBySamplesList.vcf_subset[3]

    File? outlier_filtered_manta_vcf_index = SubsetVcfBySamplesList.vcf_subset_index[0]
    File? outlier_filtered_wham_vcf_index = SubsetVcfBySamplesList.vcf_subset_index[1]
    File? outlier_filtered_melt_vcf_index = SubsetVcfBySamplesList.vcf_subset_index[2]
    File? outlier_filtered_depth_vcf_index = SubsetVcfBySamplesList.vcf_subset_index[3]

    File outlier_filtered_pesr_vcf = MergePesrVcfs.concat_vcf
    File outlier_filtered_pesr_vcf_index = MergePesrVcfs.concat_vcf_idx
    File filtered_batch_samples_file = FilterSampleList.filtered_samples_file
  }
}
