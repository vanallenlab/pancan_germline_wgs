# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Quality control of a joint-genotyped VCF


version 1.0


import "Utilities.wdl" as utils


workflow VcfQc {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    Int records_per_shard = 1000

    String output_prefix

    String bcftools_docker
    String g2c_vcf_qc_docker
  }

  # Get list of samples present in input VCFs
  call utils.GetSamplesFromVcfHeader as GetSamplesInVcf {
    input:
      vcf = vcfs[0],
      vcf_idx = vcf_idxs[0],
      bcftools_docker = bcftools_docker
  }

  # Check the size of each input vcf and shard even more if necessary
  scatter ( vcf_info in zip(vcfs, vcf_idxs) ) {

    File vcf = vcf_info.left
    File vcf_idx = vcf_info.right
    
    call utils.CountRecordsInVcf as CheckVcfSize {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        bcftools_docker = bcftools_docker
    }

    if ( CheckVcfSize.n_records > records_per_shard ) {
      call utils.ShardVcf {
        input:
          vcf = vcf,
          vcf_idx = vcf_idx,
          records_per_shard = records_per_shard,
          bcftools_docker = bcftools_docker
      }
    }

    Array[File] sharded_vcf = select_first([ShardVcf.vcf_shards, [vcf]])
    Array[File] sharded_vcf_idx = select_first([ShardVcf.vcf_shard_idxs, [vcf_idx]])
  }

  # Flatted all VCF shards into a single array for parallelization
  Array[File] vcf_shards = flatten(sharded_vcf)
  Array[File] vcf_shard_idxs = flatten(sharded_vcf_idx)

  # Scatter over all input VCF shards and collect necessary metrics for each
  scatter ( shard_info in zip(vcf_shards, vcf_shard_idxs) ) {

    File vcf_shard = shard_info.left
    File vcf_shard_idx = shard_info.right

    # TODO: add metric collection tasks here (TBD)
  }

  # Collapse all metrics across all shards
  # TODO: implement this

  # Visualize QC metrics
  # TODO: implement this

  # Package all QC plots and summary data in a single tarball
  # TODO: implement this

  output {}
}
