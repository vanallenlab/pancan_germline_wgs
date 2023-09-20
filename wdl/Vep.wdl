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

    File reference_fasta
    File vep_cache_tarball # VEP cache tarball downloaded from Ensembl
    Array[File?] other_vep_files # All other files needed for VEP. These will be moved to execution directory.
    Array[String] vep_options = [""]
    String vep_assembly = "GRCh38"
    Int vep_version = 110

    Int records_per_shard = 50000
    Boolean combine_output_vcfs = false
    String? cohort_prefix

    Boolean reindex_on_task = false # Debugging parameter

    String bcftools_docker
    String vep_docker = "vanallenlab/g2c-vep:latest"
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
      call RunVep {
        input:
          vcf = shard_info.left,
          vcf_idx = shard_info.right,
          reference_fasta = reference_fasta,
          vep_cache_tarball = vep_cache_tarball,
          other_vep_files = other_vep_files,
          vep_options = vep_options,
          vep_assembly = vep_assembly,
          vep_version = vep_version,
          reindex = reindex_on_task,
          docker = vep_docker
      }
    }
    

    call Tasks.ConcatVcfs as ConcatInnerShards {
      input:
        vcfs = RunVep.annotated_vcf,
        vcf_idxs = RunVep.annotated_vcf_idx,
        out_prefix = basename(vcf, ".vcf.gz") + ".vep",
        bcftools_docker = bcftools_docker
    }
  }

  if ( combine_output_vcfs ) {
    call Tasks.ConcatVcfs as ConcatOuterShards {
      input:
        vcfs = ConcatInnerShards.merged_vcf,
        vcf_idxs = ConcatInnerShards.merged_vcf_idx,
        out_prefix = select_first([cohort_prefix, "all_input_vcfs"]) + ".vep.merged",
        bcftools_docker = bcftools_docker
    }
  }

  output {
    Array[File] annotated_vcfs = ConcatInnerShards.merged_vcf
    Array[File] annotated_vcf_idxs = ConcatInnerShards.merged_vcf_idx
    File? combined_annotated_vcf = ConcatOuterShards.merged_vcf
    File? combined_annotated_vcf_idx = ConcatOuterShards.merged_vcf_idx
  }
}


task RunVep {
  input {
    File vcf
    File vcf_idx

    File reference_fasta
    File vep_cache_tarball
    Array[File?] other_vep_files
    String vep_assembly
    Array[String] vep_options
    Int vep_max_sv_size = 50
    Int vep_version = 110
    Boolean reindex = false

    Float mem_gb = 7.5
    Int n_cpu = 4
    Int? disk_gb

    String docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + ".vep.vcf.gz"
  Int default_disk_gb = ceil(10 * size([vcf, vep_cache_tarball, reference_fasta], "GB")) + 50

  command <<<
    set -eu -o pipefail

    # Unpack contents of cache into $VEP_CACHE/
    # Note that $VEP_CACHE is a default ENV variable set in VEP docker
    tar -xzvf ~{vep_cache_tarball} -C $VEP_CACHE/

    # Relocate other_vep_files to execution directory
    if [ ~{defined(other_vep_files)} == "true" ]; then
      while read file; do
        mv $file ./
      done < ~{write_lines(select_all(other_vep_files))}
    fi

    # Recompress & reindex VCF if optioned
    if [ ~{reindex} == "true" ]; then
      zcat ~{vcf} | bgzip -c > input.vcf.bgz
      tabix -p vcf -f input.vcf.bgz
    else
      mv ~{vcf} input.vcf.bgz
      mv ~{vcf_idx} input.vcf.bgz.tbi
    fi

    vep \
      --input_file input.vcf.bgz \
      --format vcf \
      --output_file ~{out_filename} \
      --vcf \
      --verbose \
      --force_overwrite \
      --species homo_sapiens \
      --assembly ~{vep_assembly} \
      --max_sv_size ~{vep_max_sv_size} \
      --offline \
      --no_stats \
      --cache \
      --dir_cache $VEP_CACHE/ \
      --cache_version ~{vep_version} \
      --dir_plugins $VEP_PLUGINS/ \
      --fasta ~{reference_fasta} \
      --minimal \
      --nearest gene \
      --distance 10000 \
      --numbers \
      --hgvs \
      --no_escape \
      --symbol \
      --canonical \
      --domains \
      ~{sep=" " vep_options}

    tabix -f ~{out_filename}

  >>>

  output {
    File annotated_vcf = "~{out_filename}"
    File annotated_vcf_idx = "~{out_filename}.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: 25
    preemptible: 3
  }
}
