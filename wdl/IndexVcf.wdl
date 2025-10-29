# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Index a VCF with tabix


version 1.0


import "Utilities.wdl" as Utilities
import "IndexBam.wdl" as IndexBam


workflow IndexVcf {
  input {
    File vcf
    Boolean copy_index_to_vcf_bucket = false
    String docker = "vanallenlab/g2c_pipeline:latest"
  }

  call Utilities.MakeTabixIndex as IndexVcf {
     input:
       input_file = vcf,
       docker = docker
   }

  if (copy_index_to_vcf_bucket) {
    call IndexBam.CopyIndex {
      input:
        bai = IndexVcf.tbi,
        bam = vcf,
        suffix = "tbi"
    }
  }
  
  output {
    File vcf_index = select_first([CopyIndex.bai_copy, IndexVcf.tbi])
  }
}

