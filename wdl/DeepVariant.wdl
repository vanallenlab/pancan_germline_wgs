# WDL to run DeepVariant on a single CRAM to generate single-sample gVCF
# Contact: Ryan_Collins@dfci.harvard.edu

# Based on Google's WGS case study:
# https://github.com/google/deepvariant/blob/r1.5/docs/deepvariant-case-study.md


version 1.0


workflow DeepVariant {
  input {
    File bam
    File bam_idx
    File ref_fa
    File ref_fai
    String deepvariant_docker
    String dv_model_type = "WGS"

    # DV gVCF optimization parameter to reduce output file size
    # See: https://github.com/google/deepvariant/blob/r0.5/docs/deepvariant-gvcf-support.md
    Int gq_binsize = 3

    Int? disk_size
  }

  String sample_name = basename(bam, ".bam")
    
  call RunDV { 
    input:
      bam = bam,
      bam_idx = bam_idx,
      ref_fa = ref_fa,
      ref_fai = ref_fai,
      sample_name = sample_name,
      gq_binsize = gq_binsize,
      dv_model_type = dv_model_type,
      disk_size = disk_size,
      deepvariant_docker = deepvariant_docker
  }

  output {
      File vcf = RunDV.vcf
      File vcf_idx = RunDV.vcf_idx
      File gvcf = RunDV.gvcf
  }
}


# DeepVariant one-shot invocation
# See https://github.com/google/deepvariant/blob/r1.5/docs/deepvariant-case-study.md
task RunDV {
  input {
    File bam
    File bam_idx
    File ref_fa
    File ref_fai
    String sample_name
    String dv_model_type
    String deepvariant_docker

    Int gq_binsize
    Int? disk_size
  }

  Int default_disk_size = ceil( ( 2 * size([bam, ref_fa], "GB") ) + 20.0 )

  command <<<
    set -eu -o pipefail

    mkdir -p output
    mkdir -p output/intermediate_results_dir

    /opt/deepvariant/bin/run_deepvariant \
      --model_type ~{dv_model_type} \
      --ref ~{ref_fa} \
      --reads ~{bam} \
      --output_vcf output/~{sample_name}.output.vcf.gz \
      --output_gvcf output/~{sample_name}.output.gvcf.gz \
      --num_shards $( nproc ) \
      --make_examples_extra_args gvcf_gq_binsize=~{gq_binsize} \
      --intermediate_results_dir output/intermediate_results_dir

    tabix -p vcf -f output/~{sample_name}.output.vcf.gz
  >>>

  runtime {
    docker: "~{deepvariant_docker}"
    cpu: "32"
    preemptible: 3
    disks: "local-disk " + select_first([disk_size, default_disk_size]) + " HDD"
  }

  output {
    File vcf = "output/~{sample_name}.output.vcf.gz"
    File vcf_idx = "output/~{sample_name}.output.vcf.gz.tbi"
    File gvcf = "output/~{sample_name}.output.gvcf.gz"
  }
}
