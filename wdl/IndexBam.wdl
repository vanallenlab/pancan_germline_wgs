# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Index a BAM with samtools


version 1.0


workflow IndexBam {
  input {
    File bam
    Boolean copy_index_to_bam_bucket = false
    String suffix = "bai"
    String docker = "vanallenlab/g2c_pipeline:latest"
  }

  call MakeIndex {
    input:
      bam = bam,
      docker = docker
  }

  if (copy_index_to_bam_bucket) {
    call CopyIndex {
      input:
        bai = MakeIndex.bai,
        bam = bam,
        suffix = suffix
    }
  }
  
  output {
    File bam_index = select_first([CopyIndex.bai_copy, MakeIndex.bai])
  }
}


task MakeIndex {
  input {
    File bam
    Int n_cpu = 4
    Float mem_gb = 3.75
    String docker
  }

  Int disk_size = ceil( ( 1.25 * size(bam, "GB") ) + 10 )
  String bai_fname = basename(bam) + ".bai"

  command <<<
    set -eu -o pipefail

    samtools index \
      -@ ~{n_cpu} \
      ~{bam} \
      ~{bai_fname}
  >>>

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
 }

  output {
    File bai = "~{bai_fname}"
  }
}


task CopyIndex {
  input {
    File bai
    String bam
    String suffix = "bai"
  }

  command <<<
    set -eu -o pipefail

    gsutil -m cp ~{bai} ~{bam}.~{suffix}
  >>>

  runtime {
    docker: "google/cloud-sdk:525.0.0-20250603"
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 10 HDD"
    preemptible: 3
 }

  output {
    String bai_copy = "~{bam}.~{suffix}"
  }
}

