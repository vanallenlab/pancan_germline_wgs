# The Genomic Architecture of Human Cancers
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Download fastqs for a single sample based on SRA accession (SRR ID)


version 1.0


workflow FastqFromSra {
  input {
    String accession
    File ngc
    String docker
    Int? diskGb
  }

  call DumpFastqs {
    input:
      accession = accession,
      ngc = ngc,
      docker = docker,
      diskGb = diskGb
  }

  output {
    File fastq_R1 = DumpFastqs.R1
    File fastq_R2 = DumpFastqs.R2
  }
}


task DumpFastqs {
  input {
    String accession
    File ngc
    String docker
    Int diskGb = 100
  }

  command <<<

    set -eu -o pipefail

    # Prefetch
    prefetch \
      --ngc ~{ngc} \
      --type sra \
      --max-size ~{ceil(2 * diskGb / 3)}g \
      --progress \
      ~{accession}

    # Dump fastqs
    fasterq-dump \
      --ngc ~{ngc} \
      --split-3 \
      --progress \
      --disk-limit ~{ceil(2 * diskGb / 3)}g \
      --disk-limit-tmp ~{ceil(2 * diskGb / 3)}g \
      ~{accession}

    # Compress
    cat ~{accession}_1.fastq | gzip -c > ~{accession}.R1.fq.gz
    cat ~{accession}_2.fastq | gzip -c > ~{accession}.R2.fq.gz

  >>>

  output {
    File R1 = "~{accession}.R1.fq.gz"
    File R2 = "~{accession}.R2.fq.gz"
  }

  runtime {
    cpu: 4
    memory: "7.75 GiB"
    disks: "local-disk " + diskGb + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 2
  }
}

