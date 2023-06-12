# The Genomic Architecture of Human Cancers
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Download aligned BAM for a single sample based on SRA accession (SRR ID)


version 1.0


workflow BamFromSra {
  input {
    String accession
    File? ngc
    String docker
    String? outfile_prefix
    Int? diskGb
    Int? compression
  }

  call DumpBam {
    input:
      accession = accession,
      ngc = ngc,
      docker = docker,
      outfile_prefix = outfile_prefix,
      diskGb = diskGb,
      compression = compression
  }

  output {
    File bam = DumpBam.bam_out
    File bai = DumpBam.bai_out
  }
}


task DumpBam {
  input {
    String accession
    File? ngc
    String docker
    String? outfile_prefix
    Int diskGb = 250
    Int compression = 2
  }
  String outfile_name = select_first([outfile_prefix, accession]) + ".bam"

  command <<<

    set -eu -o pipefail

    # Prefetch
    prefetch \
      ~{if defined(ngc) then "--ngc ~{ngc}" else ""} \
      --type sra \
      --max-size ~{ceil(2 * diskGb / 3)}g \
      --progress \
      ~{accession}

    # Dump SAM and save as BAM
    sam dump \
      --unaligned \
      --header \
      --verbose \
      ~{accession} \
    | samtools view \
      --bam \
      --output-fmt-option level=~{compression} \
      --write-index \
      -o ~{outfile_name}##idx##~{outfile_name}.bai

  >>>

  output {
    File bam_out = "~{accession}.bam"
    File bai_out = "~{accession}.bam.bai"
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

