# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Reblock a GVCF according to GATK-4 "biggest practices"
# See: https://broadinstitute.github.io/warp/blog/Nov21_ReblockedGVCF
# And: https://raw.githubusercontent.com/broadinstitute/warp/develop/tasks/broad/GermlineVariantDiscovery.wdl
# And: https://sites.google.com/a/broadinstitute.org/legacy-gatk-forum-discussions/2019-06-28-2019-02-11/23476-Processing-a-large-number-of-gVCFs-files-in-a-local-cluster


version 1.0


workflow ReblockGvcf {
  input {
    File input_gvcf
    File input_gvcf_index
    File ref_fasta
    File ref_fai
    File ref_dict
    File sample_name
    String gatk_docker = "us.gcr.io/broad-gatk/gatk:4.4.0.0"
  }

  call GnarlyReblock {
    input:
      gvcf = input_gvcf,
      gvcf_idx = input_gvcf_index,
      ref_fasta = ref_fasta,
      ref_fai = ref_fai,
      ref_dict = ref_dict,
      out_fname = sample_name + ".reblocked.g.vcf.gz",
      docker = gatk_docker
  }
  
  output {
    File reblocked_gvcf = GnarlyReblock.output_gvcf
    File reblocked_gvcf_idx = GnarlyReblock.output_gvcf_idx
  }
}


task GnarlyReblock {
  input {
    File gvcf
    File gvcf_idx
    File ref_fasta
    File ref_fai
    File ref_dict
    String out_fname
    String docker
  }

  Int disk_size = ceil((size(gvcf, "GiB")) * 4) + 20

  command <<<
    set -eu -o pipefail

    gatk --java-options "-Xms3000m -Xmx3000m" \
      ReblockGVCF \
      -R ~{ref_fasta} \
      -V ~{gvcf} \
      -drop-low-quals \
      -rgq-threshold 10 \
      -do-qual-approx \
      --floor-blocks -GQB 20 -GQB 30 -GQB 40 \
      --create-output-variant-index \
      -O ~{out_fname}
  >>>

  runtime {
    memory: "3750 MiB"
    disks: "local-disk " + disk_size + " HDD"
    bootDiskSizeGb: 15
    preemptible: 3
    docker: docker
  }

  output {
    File output_gvcf = out_fname
    File output_gvcf_idx = out_fname + ".tbi"
  }
}

