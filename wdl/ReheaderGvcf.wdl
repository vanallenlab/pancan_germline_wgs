# The Genomic Architecture of Human Cancers
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Reheader and rename a gVCF


version 1.0


workflow ReheaderGvcf {
  input {
    File input_gvcf
    File sample_name
    String bcftools_docker
    Int? bgzip_compression
  }

  call Reheader {
    input:
      input_gvcf = input_gvcf,
      sample_name = sample_name,
      bgzip_compression = bgzip_compression,
      docker = bcftools_docker
  }
  
  output {
    File output_gvcf = Reheader.gvcf_out
    File output_gvcf_idx = Reheader.gvcf_out_idx
  }
}


task Reheader {
  input {
    File input_gvcf
    String sample_name
    Int bgzip_compression = -1
    String docker
  }
  Int disk_size = ceil( ( 2.3 * size(input_gvcf, "GB") ) + 10.0 )

  command <<<
    set -eu -o pipefail

    bcftools reheader \
      --samples ~{write_lines([sample_name])} \
      ~{input_gvcf} \
    | gunzip -c \
    | bgzip -l ~{bgzip_compression} -c \
    > ~{sample_name}.g.vcf.gz
    tabix -p vcf -f ~{sample_name}.g.vcf.gz
  >>>

  runtime {
    docker: docker
    memory: "3.5 GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
 }

  output {
    File gvcf_out = "~{sample_name}.g.vcf.gz"
    File gvcf_out_idx = "~{sample_name}.g.vcf.gz.tbi"
  }
}
