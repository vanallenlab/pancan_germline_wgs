# The Genomic Architecture of Human Cancers
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Convert a paired-end BAM to paired, gzipped fastqs


version 1.0


workflow Bam2Fastq {
  input {
    File bam
    File sample_name
    String samtools_docker
    Int? cpu
  }

  call B2F {
    input:
      bam = bam,
      sample_name = sample_name,
      docker = samtools_docker,
      cpu = cpu
  }
  
  output {
    File fastq_1 = B2F.fq1
    File fastq_2 = B2F.fq2
  }
}


task B2F {
  input {
    File bam
    String sample_name
    String docker
    Int cpu = 4
    Float mem_gb_per_cpu = 1.75
    Int? disk_gb
  }

  Int default_disk_size = ceil( ( 4 * size(bam, "GB") ) + 10.0 )
  Float total_mem_gb = mem_gb_per_cpu * cpu

  command <<<
    set -eu -o pipefail

    samtools sort -n -@ ~{cpu} ~{bam} \
    | samtools fastq -c 7 \
      -1 ~{sample_name}.R1.fq.gz \
      -2 ~{sample_name}.R2.fq.gz \
      -@ ~{cpu}
  >>>

  runtime {
    docker: docker
    memory: "~{total_mem_gb} GB"
    disks: "local-disk " + select_first([disk_gb, default_disk_size]) + " HDD"
    preemptible: 3
  }

  output {
    File fq1 = "~{sample_name}.R1.fq.gz"
    File fq2 = "~{sample_name}.R2.fq.gz"
  }
}
