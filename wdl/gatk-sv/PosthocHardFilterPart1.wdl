# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Post hoc SV site-level hard filters for raw GATK-SV output (from module 15/16)

# Part 1 of 2, as follows:
# 1. Basic site-level hard filters
# [Manual step between] Outlier identification & exclusion
# 2. Outlier sample exclusion & second layer of GQ-dependent hard filters


version 1.0


import "Utilities.wdl" as utils


workflow PosthocHardFilterPart1 {
  input {
    File vcf
    File vcf_idx
    String bcftools_docker
  }

  call HardFilterPart1 {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      docker = bcftools_docker
  }

  output {
    File filtered_vcf = HardFilterPart1.filtered_vcf
    File filtered_vcf_idx = HardFilterPart1.filtered_vcf_idx
  }
}


task HardFilterPart1 {
  input {
    File vcf
    File vcf_idx

    String docker
    Float mem_gb = 7.5
    Int n_cpu = 4
  }

  Int disk_gb = ceil(2 * size([vcf], "GB")) + 10
  String outfile = basename(vcf, ".vcf.gz") + ".posthoc_filtered.vcf.gz"

  command <<<
    set -euo pipefail

    # Post hoc filters as follows:
    # 1. Remove all unresolved variants (breakends)
    # 2. Remove all deletions called only by the Wham algorithm (>98% FDR)

    bcftools view \
      --include '(FILTER="PASS" | FILTER="MULTIALLELIC")' \
      -Oz -o "~{outfile}" \
      ~{vcf}

    tabix -p vcf "~{outfile}"
  >>>

  output {
    File filtered_vcf = "~{outfile}"
    File filtered_vcf_idx = "~{outfile}.tbi"
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

