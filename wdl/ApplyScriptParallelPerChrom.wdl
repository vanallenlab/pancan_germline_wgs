# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Generic, flexible WDL to scatter a user-supplied script over all chromosomes in one or more VCF(s)


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/Utilities.wdl" as Utilities


workflow ApplyScriptParallelPerChrom {
  input {
    File vcf
    File script                      # User-supplied script staged in a google bucket. Must take input and output VCFs as final two positional arguments 
    String? out_vcf_prefix

    String exec_prefix = "source"    # Command to interpret script as an executable
    String script_options = ""       # Any other command-line options to be passed to the script
    Array[String?] script_files      # Files needed by script to be localized to execution directory

    File? vcf_idx                    # Recommended but not strictly required
    File? ref_fai                    # Used to determine contigs for parallelization. By default will take all contigs in vcf header
    String? bcftools_concat_options 

    Float? script_mem_gb
    Int? script_cpu_cores
    Int? script_disk_gb
    Float? merge_mem_gb
    Int? merge_cpu_cores
    Int? merge_disk_gb

    String bcftools_docker            # Any linux-based image with bcftools & tabix installed. Must also have dependencies for running user-supplied script.
  }

  # Get contigs for scatter depending on user input
  if (defined(ref_fai)) {
    call Utilities.GetContigsFromFai {
      input:
        ref_fai = select_first([ref_fai]),
        docker = bcftools_docker
    }
  }
  if (!defined(ref_fai)) {
    call Utilities.GetContigsFromVcfHeader {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        docker = bcftools_docker
    }
  }
  Array[String] contigs = select_first([GetContigsFromFai.contigs, 
                                        GetContigsFromVcfHeader.contigs])

  # Apply script to each chromosome in parallel
  scatter (contig in contigs) {
    call ApplyScript {
      input:
        vcf = vcf,
        vcf_idx = select_first([vcf_idx, GetContigsFromVcfHeader.vcf_idx_out]),
        contig = contig,
        script = script,
        exec_prefix = exec_prefix,
        script_options = script_options,
        script_files = script_files,
        mem_gb = script_mem_gb,
        cpu_cores = script_cpu_cores,
        disk_gb = script_disk_gb,
        docker = bcftools_docker

    }
  }

  # Merge outputs from per-chromosome scatter
  call Utilities.ConcatVcfs {
    input:
      vcfs = ApplyScript.vcf_out,
      vcf_idxs = ApplyScript.vcf_out_idx,
      out_prefix = select_first([out_vcf_prefix, basename(vcf, ".vcf.gz") + ".processed"]),
      bcftools_concat_options = bcftools_concat_options,
      mem_gb = merge_mem_gb,
      cpu_cores = merge_cpu_cores,
      disk_gb = merge_disk_gb,
      bcftools_docker = bcftools_docker
  }

  output {
    File output_vcf = ConcatVcfs.merged_vcf
    File output_vcf_idx = ConcatVcfs.merged_vcf_idx
  }
}


task ApplyScript {
  input {
    File vcf
    File? vcf_idx
    String? contig
    File script
    String exec_prefix
    String script_options
    Array[String?] script_files

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String docker
  }

  String vcf_basename = basename(vcf, ".vcf.gz")
  String fname_tag = if defined(contig) then "." + contig else ""
  String contig_str = if defined(contig) then contig else ""
  String in_vcf_name = vcf_basename + fname_tag + ".in.vcf.gz"
  String out_vcf_name = vcf_basename + fname_tag + ".out.vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcf, "GB")) + 10
  Int use_disk_gb = select_first([disk_gb, default_disk_gb])

  command <<<
    set -eu -o pipefail

    if [ ~{length(script_files)} -gt 0 ]; then
      cat ~{write_lines(select_all(script_files))} \
      | gsutil -m cp -I ./
    fi

    if [ ~{defined(vcf_idx)} == "false" ]; then
      tabix -p vcf -f ~{vcf}
    fi

    # Subset VCF to contig of interest, if optioned
    if [ ~{defined(contig)} == "true" ]; then
      tabix -h ~{vcf} "~{contig_str}" | bgzip -c > "~{in_vcf_name}"
      tabix -p vcf -f "~{in_vcf_name}"
    else
      mv ~{vcf} ~{in_vcf_name}
      mv ~{vcf}.tbi ~{in_vcf_name}.tbi
    fi

    # Execute script
    cmd="~{exec_prefix} ~{script} ~{script_options} ~{in_vcf_name} ~{out_vcf_name}"
    eval "$cmd"
    tabix -p vcf -f "~{out_vcf_name}"
  >>>

  output {
    File vcf_out = "~{out_vcf_name}"
    File vcf_out_idx = "~{out_vcf_name}.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: cpu_cores
    disks: "local-disk " + use_disk_gb + " HDD"
    preemptible: 3
  }
}
