# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Post hoc SV site-level hard filters for raw GATK-SV output (from module 15/16)


version 1.0


import "Utilities.wdl" as utils


workflow PosthocHardFilter {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String bcftools_docker
  }

  Array[Pair[File, File]] vcf_infos = zip(vcfs, vcf_idxs)

  scatter ( vcf_info in vcf_infos ) {

    File vcf = vcf_info.left
    File vcf_idx = vcf_info.right

    call HardFilter {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        docker = bcftools_docker
    }
  }

  output {
    Array[File] filtered_vcfs = HardFilter.filtered_vcf
    Array[File] filtered_vcf_idxs = HardFilter.filtered_vcf_idx
  }
}


task HardFilter {
  input {
    File vcf
    File vcf_idx

    String docker
    Float mem_gb = 3.75
    Int n_cpu = 2
  }

  Int disk_gb = ceil(2 * size([vcf], "GB")) + 10
  String outfile = basename(vcf, ".vcf.gz") + ".posthoc_filtered.vcf.gz"

  command <<<
    set -euo pipefail

    # Post hoc filters as follows:
    # 1. Remove all unresolved variants (breakends)
    # 2. Remove all deletions called only by the Wham algorithm (>98% FDR)
    # 3. Exclude variants with no non-reference genotypes with GQ > 1

    # First, find all variant IDs with at least one sample with non-ref GQ > 1
    bcftools query -i '(GT="alt" & FORMAT/GQ>1) | FILTER="MULTIALLELIC"' -f '%ID\n' ~{vcf} \
    | sort -V | uniq > pass_vids.list

    # Second, restrict to passing variants above and apply all other filters
    bcftools view \
      --include 'ID=@pass_vids.list & (FILTER="PASS" | FILTER="MULTIALLELIC")' \
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

