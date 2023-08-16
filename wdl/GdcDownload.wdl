# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Download a single file from NCI GDC using gdc-client
# Does not perform any additional processing beyond strictly downloading the file


version 1.0


workflow GdcDownload {
  input {
    String uuid
    String filename
    File? token
    String docker
    Int? diskGb
  }

  call DownloadFile {
    input:
      uuid = uuid,
      filename = filename,
      token = token,
      docker = docker,
      diskGb = diskGb
  }

  output {
    File outfile = DownloadFile.file_out
  }
}


task DownloadFile {
  input {
    String uuid
    String filename
    File? token
    String docker
    Int diskGb = 100
  }

  command <<<

    set -eu -o pipefail

    # Download using gdc-client
    /opt/gdc-client download \
      ~{if defined(token) then "--token " + select_first([token]) else ""} \
      ~{uuid}

  >>>

  output {
    File file_out = "~{filename}"
  }

  runtime {
    cpu: 2
    memory: "3.75 GiB"
    disks: "local-disk " + diskGb + " HDD"
    bootDiskSizeGb: 10
    docker: docker
    preemptible: 3
    maxRetries: 2
  }
}

