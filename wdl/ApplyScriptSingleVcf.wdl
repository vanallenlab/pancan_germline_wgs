# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Apply a user-supplied script to a single input VCF with no sharding

# For an equivalent version of this script parallelized by chromosome, see:
# ApplyScriptParallelPerChrom.wdl

# For an equivalent version of this applied over multiple VCFs (with or without chromosome-parallelization), see:
# MultiVcfParallelizeScript.wdl


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/ApplyScriptParallelPerChrom.wdl" as Parallelize
import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/Utilities.wdl" as Utilities


workflow ApplyScriptSingleVcf {
  input {
    File vcf
    File? vcf_idx

    File script
    String exec_prefix = "source"
    String script_options = ""
    Array[String?] script_files

    String? out_vcf_prefix

    Float? script_mem_gb
    Int? script_cpu_cores
    Int? script_disk_gb

    String bcftools_docker
  }

  # Index input VCF if necessary
  if (!defined(vcf_idx)) {
    call Utilities.IndexVcf {
      input:
        vcf = vcf,
        docker = bcftools_docker
    }
  }

  call Parallelize.ApplyScript as ApplyScript {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      script = script,
      exec_prefix = exec_prefix,
      script_options = script_options,
      script_files = script_files,
      mem_gb = script_mem_gb,
      cpu_cores = script_cpu_cores,
      disk_gb = script_disk_gb,
      docker = bcftools_docker
  }

  output {
    File output_vcf = ApplyScript.vcf_out
    File output_vcf_idx = ApplyScript.vcf_out_idx
  }
}

