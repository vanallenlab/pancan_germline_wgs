# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# WDL to run VEP on one or more input VCFs


version 1.0


workflow Vep {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    File reference_fasta
    File vep_cache_tarball # VEP cache tarball downloaded from Ensembl
    Array[String?] vep_remote_files # URIs for files stored in Google buckets that can be remotely sliced using tabix for each shard
    Array[File?] vep_remote_file_indexes # Indexes corresponding to vep_remote_files
    Array[File?] other_vep_files # All other files needed for VEP. These will be localized in full to each VM and moved to execution directory.

    Array[String] vep_options = [""]
    String vep_assembly = "GRCh37"
    Int vep_version = 110

    Boolean shard_vcfs = false
    Int remote_query_buffer = 2
    String? cohort_prefix

    String bcftools_docker
    String vep_docker = "vanallenlab/g2c-vep:latest"
  }

  Int n_remote_files = length(select_all(vep_remote_files))
  Boolean any_remote = n_remote_files > 0



  Array[File?] all_other_vep_files = flatten(select_all([other_vep_files]))

  call RunVep {
    input:
      vcf = vcfs[0],
      vcf_idx = vcf_idxs[0],
      reference_fasta = reference_fasta,
      vep_cache_tarball = vep_cache_tarball,
      other_vep_files = all_other_vep_files,
      vep_options = vep_options,
      vep_assembly = vep_assembly,
      vep_version = vep_version,
      docker = vep_docker
  }


  output {
    File annotated_vcfs = RunVep.annotated_vcf
    File annotated_vcf_idxs = RunVep.annotated_vcf_idx
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
    #if [ ~{defined(other_vep_files)} == "true" ]; then
    #  while read file; do
    #    mv "$file" ./
    #  done < ~{write_lines(select_all(other_vep_files))}
    #fi


    vep \
      --input_file ~{vcf} \
      --format vcf \
      --output_file ~{out_filename} \
      --vcf \
      --verbose \
      --compress_output bgzip \
      --force_overwrite \
      --species homo_sapiens \
      --assembly ~{vep_assembly} \
      --max_sv_size ~{vep_max_sv_size} \
      --offline \
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