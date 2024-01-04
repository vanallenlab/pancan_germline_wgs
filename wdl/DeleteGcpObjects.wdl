# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Parallel deletion of objects stored in GCP buckets


version 1.0


workflow DeleteGcpObjects {
  input {
    File uri_list
    Int n_cpu = 4
  }

  call DeleteObjects {
    input:
      uri_list = uri_list,
      n_cpu = n_cpu
  }
}


task DeleteObjects {
  input {
    File uri_list
    Int n_cpu
  }

  Int n_threads = 2 * n_cpu
  Float mem_gb = 2 + ( 0.5 * n_cpu )
  Int disk_gb = 10 + ceil(2 * size(uri_list, "GB"))

  command <<<
    set -eu -o pipefail

    cat ~{uri_list} \
    | gsutil \
        -m \
        -o GSUtil:parallel_process_count=1 \
        -o GSUtil:parallel_thread_count=~{n_threads} \
        rm -fI
  >>>

  runtime {
    docker: "google/cloud-sdk:latest"
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 5
 }
}

