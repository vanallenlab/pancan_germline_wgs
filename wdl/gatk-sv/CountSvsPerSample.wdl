# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Count SVs per sample per SV type
# Supports chromosome-scattered inputs for parallelization


version 1.0


workflow CountSvsPerSample {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String output_prefix
    String sv_pipeline_docker
    String g2c_pipeline_docker
  }

  Array[Pair[File, File]] vcf_infos = zip(vcfs, vcf_idxs)

  scatter ( vcf_info in vcf_infos ) {
    
    File vcf = vcf_info.left
    File vcf_idx = vcf_info.right

    call CountSvs {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        docker = sv_pipeline_docker
    }
  }

  if ( length(vcfs) > 1 ) {
    call SumCounts {
      input:
        count_tsvs = CountSvs.counts_tsv,
        output_prefix = output_prefix,
        docker = g2c_pipeline_docker
    }
  }

  output {
    File sv_counts = select_first([SumCounts.summed_tsv, CountSvs.counts_tsv[0]])
  }
}


task CountSvs {
  input {
    File vcf
    File vcf_idx

    String docker
    Float mem_gb = 3.75
    Int n_cpu = 2
  }

  Int disk_gb = ceil(2 * size([vcf], "GB")) + 10
  String outfile = basename(vcf, ".vcf.gz") + ".counts.tsv"

  command <<<
    set -euo pipefail

    svtk count-svtypes ~{vcf} ~{outfile}
  >>>

  output {
    File counts_tsv = "~{outfile}"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}


task SumCounts {
  input {
    Array[File] count_tsvs
    String output_prefix

    String docker
    Float mem_gb = 3.75
    Int n_cpu = 2
  }

  Int disk_gb = ceil(2 * size(count_tsvs, "GB")) + 10
  String outfile = output_prefix + ".counts.tsv"

  command <<<
    set -euo pipefail

    # TBD: python script packaged in g2c docker here
  >>>

  output {
    File summed_tsv = "~{outfile}"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}
