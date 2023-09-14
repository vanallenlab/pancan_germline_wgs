# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# WDL to run VEP on one or more input VCFs


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/Utilities.wdl" as Tasks


workflow Vep {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    Int records_per_shard = 5000
    Boolean combine_output_vcfs = false
    String? cohort_prefix

    String bcftools_docker
  }

  scatter ( vcf_info in zip(vcfs, vcf_idxs) ) {
    File vcf = vcf_info.left
    File vcf_idx = vcf_info.right

    call Tasks.ShardVcf {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        records_per_shard = records_per_shard,
        bcftools_docker = bcftools_docker
    }

    scatter ( shard_info in zip(ShardVcf.vcf_shards, ShardVcf.vcf_shard_idxs) ) {
      # TODO: call VEP here
    }
    

    call Tasks.ConcatVcfs as ConcatInnerShards {
      input:
        vcfs = RunVep.annotated_vcf,
        vcf_idxs = RunVep.annotated_vcf_idx,
        out_prefix = basename(vcf, ".vcf.gz") + ".vep",
        docker = bcftools_docker
    }
  }

  if ( combine_output_vcfs ) {
    call ConcatVcfs as ConcatOuterShards {
      input:
        vcfs = ConcatInnerShards.merged_vcf,
        vcf_idxs = ConcatInnerShards.merged_vcf_idx,
        out_prefix = select_first([cohort_prefix, "all_input_vcfs"]) + ".vep.merged"),
        docker = bcftools_docker
    }
  }

  output {
    Array[File] annotated_vcfs = ConcatInnerShards.merged_vcf
    Array[File] annotated_vcf_idxs = ConcatInnerShards.merged_vcf_idx
    File? combined_annotated_vcf = ConcatOuterShards.merged_vcf
    File? combined_annotated_vcf_idx = ConcatOuterShards.merged_vcf_idx
  }
}


