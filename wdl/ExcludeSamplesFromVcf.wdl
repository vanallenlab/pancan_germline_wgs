# The Genomic Architecture of Human Cancers
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Simple WDL to exclude samples from a VCF using BCFTools


version 1.0


workflow ExcludeSamplesFromVcf {
  input {
    File vcf
    File exclude_samples_list
    
    String? outfile_prefix
    Int? compression_level
    String? bcftools_filter_options
    Boolean update_info = true

    String docker
  }

  call ExcludeSamples {
    input:
      vcf = vcf,
      exclude_samples_list = exclude_samples_list,
      outfile_prefix = outfile_prefix,
      compression_level = compression_level,
      bcftools_filter_options = bcftools_filter_options,
      update_info = update_info,
      docker = docker
  }

  output {
    File filtered_vcf = ExcludeSamples.out_vcf
    File filtered_vcf_idx = ExcludeSamples.out_vcf_idx
  }
}


task ExcludeSamples {
  input {
    File vcf
    File exclude_samples_list

    String? outfile_prefix
    Int compression_level = 2
    String bcftools_filter_options = ""
    Boolean update_info

    String docker
  }

  String default_prefix = basename(vcf, ".vcf.gz") + ".samples_excluded"
  String outfile = select_first([outfile_prefix, default_prefix]) + ".vcf.gz"

  Int disk_gb = ceil(2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools view \
      --compression-level ~{compression_level} \
      --output ~{outfile} \
      --output-type z \
      --threads 2 \
      --samples-file ^~{exclude_samples_list} \
      --force-samples \
      ~{if update_info then "" else "--no-update"} \
      ~{bcftools_filter_options} \
      ~{vcf}
    tabix -p vcf -f ~{outfile}
  >>>

  output {
    File out_vcf = "~{outfile}"
    File out_vcf_idx = "~{outfile}.tbi"
  }

  runtime {
    docker: docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}
