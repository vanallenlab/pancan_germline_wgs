# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Post hoc outlier sample exclusion and site-level hard filtering for GATK-HC output

# Part 1 of 2, as follows:
# 1. Variant normalization & basic site-level hard filters
# [Manual step between] Outlier identification
# 2. Outlier sample exclusion & reapplication of basic hard filters


version 1.0


workflow PosthocCleanupPart2 {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    File exclude_samples_list
    String bcftools_docker
  }

  Array[Pair[File, File]] vcf_infos = zip(vcfs, vcf_idxs)

  scatter ( vcf_info in vcf_infos ) {
    call CleanupPart2 {
      input:
        vcf = vcf_info.left,
        vcf_idx = vcf_info.right,
        exclude_samples_list = exclude_samples_list,
        docker = bcftools_docker
    }
  }

  output {
    Array[File] filtered_vcfs = CleanupPart2.filtered_vcf
    Array[File] filtered_vcf_idxs = CleanupPart2.filtered_vcf_idx
  }
}


task CleanupPart2 {
  input {
    File vcf
    File vcf_idx
    File exclude_samples_list

    String docker
    Float mem_gb = 3.75
    Int n_cpu = 2
  }

  Int disk_gb = ceil(2.25 * size([vcf], "GB")) + 10
  String outfile = basename(vcf, ".vcf.gz") + ".posthoc_filtered.vcf.gz"

  command <<<
    set -euo pipefail

    # Post hoc filters as follows:
    # 1. Exclude outlier samples
    # 2. Exclude variants with no non-reference genotypes with DP > 10 and GQ > 20

    bcftools view \
      --samples-file "^~{exclude_samples_list}" \
      ~{vcf} \
    | bcftools view \
      --include 'INFO/AC>0 & FORMAT/DP>10 & FORMAT/GQ>20 & GT="alt" & alt[0] != "*"' \
      -Oz -o "~{outfile}"

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

