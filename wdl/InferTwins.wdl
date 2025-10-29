# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Infer identical twins and sample replicate pairs from joint genotyped VCFs


version 1.0


import "Utilities.wdl" as Utils


workflow InferTwins {
	input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String out_prefix

    Float min_af = 0.01
    Float min_hwe_p = 0.00000005
    Float max_missing = 0.1
    Boolean only_snvs = false
    Boolean only_biallelic = true
    Float min_kin = 0.4

    String bcftools_docker
    String plink2_docker
  }

  Array[Pair[File, File]] vcf_infos = zip(vcfs, vcf_idxs)

  scatter ( vcf_info in vcf_infos ) {
    call PrepVcf {
      input:
        vcf = vcf_info.left,
        vcf_idx = vcf_info.right,
        min_af = min_af,
        min_hwe_p = min_hwe_p,
        max_missing = max_missing,
        only_snvs = only_snvs,
        only_biallelic = only_biallelic,
        docker = bcftools_docker
    }
  }

  call Utils.ConcatVcfs {
    input:
      vcfs = PrepVcf.filtered_vcf,
      vcf_idxs = PrepVcf.filtered_vcf_idx,
      out_prefix = out_prefix + ".filtered.merged",
      bcftools_docker = bcftools_docker
  }

  # Infer relatives with KING
  call InferRelatives {
    input:
      vcf = ConcatVcfs.merged_vcf,
      vcf_idx = ConcatVcfs.merged_vcf_idx,
      out_prefix = out_prefix,
      min_kin = min_kin,
      docker = plink2_docker
  }

  output {
    File king_metrics = InferRelatives.king_metrics
  }
}


task PrepVcf {
  input {
    File vcf
    File vcf_idx

    Float min_af
    Float min_hwe_p
    Float max_missing
    Boolean only_snvs
    Boolean only_biallelic

    Int? disk_gb
    Float mem_gb = 4.0
    Int n_cpu = 2
    String docker
  }

  String snv_cmd = if only_snvs then "--types snps" else ""
  String biallelic_cmd = if only_biallelic then "-m2 -M2" else ""

  String outfile = basename(vcf, ".vcf.gz") + ".filtered.vcf.gz"

  Int default_disk_gb = ceil(2.5 * size([vcf], "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools +fill-tags ~{vcf} -- -t AC,AN,AF,HWE,F_MISSING \
    | bcftools view \
      --apply-filters "PASS,." \
      --min-af ~{min_af} \
      --include "INFO/HWE >= ~{min_hwe_p} & INFO/F_MISSING <= ~{max_missing}" \
      ~{snv_cmd} \
      ~{biallelic_cmd} \
      -Oz -o ~{outfile}

    tabix -p vcf -f ~{outfile}
  >>>

  output {
    File filtered_vcf = "~{outfile}"
    File filtered_vcf_idx = "~{outfile}.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 2
    max_retries: 1
  }
}


task InferRelatives {
  input {
    File vcf
    File vcf_idx
    String out_prefix

    Float min_kin

    Float mem_gb = 32.0
    Int n_cpu = 16
    Int? disk_gb
    String docker
  }

  Int default_disk_gb = ceil((5 * size(vcf, "GB"))) + 10

  command <<<
    set -eu -o pipefail

    # Convert VCF to PLINK format
    plink2 --vcf ~{vcf} --make-bed --out ~{out_prefix}

    # Run KING
    plink2 \
      --bfile ~{out_prefix} \
      --make-king-table counts \
      --king-table-filter 0.4 \
      --out ~{out_prefix}

    # Compress KING table
    gzip -f ~{out_prefix}.kin0
  >>>

  output {
    File king_metrics = "~{out_prefix}.kin0.gz"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 2
    max_retries: 1
  }  
}