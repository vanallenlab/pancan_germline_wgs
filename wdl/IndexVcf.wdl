# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Index a VCF with tabix


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/Utilities.wdl" as Utilities
import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/IndexBam.wdl" as IndexBam


workflow IndexVcf {
  input {
    File vcf
    Boolean copy_index_to_bam_bucket = false
    String docker = "vanallenlab/g2c_pipeline:latest"
  }

  call Utilities.IndexVcf {
     input:
       vcf = vcf,
       docker = bcftools_docker
   }

  if (copy_index_to_bam_bucket) {
    call IndexBam.CopyIndex {
      input:
        bai = IndexVcf.vcf_idx,
        bam = vcf,
        suffix = "tbi",
        docker = docker
    }
  }
  
  output {
    File bam_index = select_first([CopyIndex.bai_copy, IndexVcf.vcf_idx])
  }
}

