# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Generic, flexible WDL to scatter a user-supplied script over all chromosomes in a VCF


version 1.0


workflow ApplyScriptParallelPerChrom {
  input {
    File vcf
    File script                     # User-supplied script staged in a google bucket. Must take input and output VCFs as final two positional arguments 
    String? out_vcf_prefix

    String exec_prefix = "source"   # Command to interpret script as an executable
    String script_options = ""      # Any other command-line options to be passed to the script

    File? vcf_idx                   # Recommended but not strictly required
    File? ref_fai                   # Used to determine contigs for parallelization. By default will take all contigs in vcf header
    String? bcftools_concat_options 

    Float? script_mem_gb
    Int? script_cpu_cores
    Int? script_disk_gb
    Float? merge_mem_gb
    Int? merge_cpu_cores
    Int? merge_disk_gb

    String bcftools_docker          # Any linux-based image with bcftools & tabix installed. Must also have dependencies for running user-supplied script.
  }

  # Get contigs for scatter depending on user input
  if (defined(ref_fai)) {
    call GetContigsFromFai {
      input:
        ref_fai = select_first([ref_fai]),
        docker = bcftools_docker
    }
  }
  if (!defined(ref_fai)) {
    call GetContigsFromVcfHeader {
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
    call ApplyScriptSingleContig {
      input:
        vcf = vcf,
        vcf_idx = select_first([vcf_idx, GetContigsFromVcfHeader.vcf_idx_out]),
        contig = contig,
        script = script,
        exec_prefix = exec_prefix,
        script_options = script_options,
        mem_gb = script_mem_gb,
        cpu_cores = script_cpu_cores,
        disk_gb = script_disk_gb,
        docker = bcftools_docker

    }
  }

  # Merge outputs from per-chromosome scatter
  call ConcatVcfs {
    input:
      vcfs = ApplyScriptSingleContig.vcf_out,
      vcf_idxs = ApplyScriptSingleContig.vcf_out_idx,
      out_prefix = select_first([out_vcf_prefix, basename(vcf, ".vcf.gz") + ".processed"]),
      bcftools_concat_options = bcftools_concat_options,
      mem_gb = merge_mem_gb,
      cpu_cores = merge_cpu_cores,
      disk_gb = merge_disk_gb,
      docker = bcftools_docker
  }

  output {
    File output_vcf = ConcatVcfs.merged_vcf
    File output_vcf_idx = ConcatVcfs.merged_vcf_idx
  }
}


task GetContigsFromFai {
  input {
    File ref_fai
    String docker
  }

  command <<<
    set -eu -o pipefail

    cut -f1 ~{ref_fai} > contigs.list
  >>>

  output {
    Array[String] contigs = read_lines("contigs.list")
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 10 HDD"
    preemptible: 3
  }
}


task GetContigsFromVcfHeader {
  input {
    File vcf
    File? vcf_idx
    String docker
  }

  Int disk_gb = ceil(1.2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    if [ ~{defined(vcf_idx)} == "false" ]; then
      tabix -p vcf -f ~{vcf}
    fi

    tabix -H ~{vcf} \
    | fgrep "##contig" \
    | sed 's/ID=/\t/g' \
    | cut -f2 \
    | cut -f1 -d, \
    | sort -V \
    > contigs.list
  >>>

  output {
    Array[String] contigs = read_lines("contigs.list")
    File vcf_idx_out = select_first([vcf_idx, vcf + ".tbi"])
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task ApplyScriptSingleContig {
  input {
    File vcf
    File? vcf_idx
    String contig
    File script
    String exec_prefix
    String script_options

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String docker
  }

  String vcf_basename = basename(vcf, ".vcf.gz")
  String in_vcf_name = vcf_basename + "." + contig + ".in.vcf.gz"
  String out_vcf_name = vcf_basename + "." + contig + ".out.vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    if [ ~{defined(vcf_idx)} == "false" ]; then
      tabix -p vcf -f ~{vcf}
    fi

    # Subset VCF to contig of interest
    tabix -h ~{vcf} "~{contig}" | bgzip -c > "~{in_vcf_name}"

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
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}


task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String out_prefix

    String bcftools_concat_options = ""

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String docker
  }

  String out_filename = out_prefix + ".vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcfs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools concat \
      ~{bcftools_concat_options} \
      --file-list ~{write_lines(vcfs)} \
      -O z \
      -o ~{out_filename} \
      --threads ~{cpu_cores}

    tabix -p vcf -f ~{out_filename}
  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: cpu_cores
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}