# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Common WDL tasks shared across workflows


version 1.0


task ShardVcf {
  input {
    File vcf
    File vcf_idx
    Int records_per_shard
    String bcftools_docker
  }

  String out_prefix = basename(vcf, ".vcf.gz") + ".sharded"
  Int disk_gb = ceil(3 * size(vcf, "GB"))

  command <<<
    set -eu -o pipefail

    # Make an empty shard in case the input VCF is totally empty
    bcftools view -h ~{vcf} | bgzip -c > "~{out_prefix}.0.vcf.gz"

    bcftools +scatter \
      -O z3 -o . -p "~{out_prefix}". \
      -n ~{records_per_shard} \
      ~{vcf}

    # Print all VCFs to stdout for logging purposes
    find ./ -name "*.vcf.gz"

    # Index all shards
    find ./ -name "~{out_prefix}.*.vcf.gz" \
    | xargs -I {} tabix -p vcf -f {}
  >>>

  output {
    Array[File] vcf_shards = glob("~{out_prefix}.*.vcf.gz")
    Array[File] vcf_shard_idxs = glob("~{out_prefix}.*.vcf.gz.tbi")
  }

  runtime {
    cpu: 2
    memory: "3.75 GiB"
    disks: "local-disk " + disk_gb + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: bcftools_docker
    preemptible: 3
    maxRetries: 1
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

    String bcftools_docker
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
    docker: bcftools_docker
    memory: mem_gb + " GB"
    cpu: cpu_cores
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
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


task IndexVcf {
  input {
    File vcf
    String docker
  }

  Int disk_gb = ceil(1.25 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    tabix -p vcf -f ~{vcf}
  >>>

  output {
    File vcf_idx = "~{vcf}.tbi"
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task ConcatTextFiles {
  input {
    Array[File] shards
    String concat_command = "cat"
    String? sort_command
    String? compression_command
    Boolean input_has_header = false
    String output_filename

    String docker
  }

  Int disk_gb = ceil(2 * size(shards, "GB")) + 25
  String sort = if defined(sort_command) then " | " + sort_command else ""
  String compress = if defined(compression_command) then " | " + compression_command else ""
  String posthoc_cmds = sort + " | cat header.txt - " + compress

  command <<<
    set -eu -o pipefail

    if [ "~{input_has_header}" == "true" ]; then
      ~{concat_command} ~{shards[0]} \
      | head -n1 > header.txt || true
    else
      touch header.txt
    fi

    ~{concat_command} ~{sep=" " shards} ~{posthoc_cmds} > ~{output_filename}
  >>>

  output {
    File merged_file = "~{output_filename}"
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}