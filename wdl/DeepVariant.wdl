# WDL to run DeepVariant on a single CRAM to generate single-sample gVCF
# Contact: Ryan_Collins@dfci.harvard.edu

# Based on DNANexus' original (archived) workflow. See:
# https://github.com/dnanexus-rnd/DeepVariant-GLnexus-WDL/blob/master/wdl/DeepVariant.wdl


version 1.0


workflow DeepVariant {
  input {
    File bam
    File bam_idx
    File ref_fa
    File ref_fai
    String deepvariant_docker
    String deepvariant_model_path = "/opt/models/wgs/model.ckpt"

    # DV gVCF optimization parameter to reduce output file size
    # See: https://github.com/google/deepvariant/blob/r0.5/docs/deepvariant-gvcf-support.md
    Int? gq_binsize = 3
  }

  String sample_name = basename(bam, ".bam")
    
  call make_examples { 
    input:
      bam = bam,
      bam_idx = bam_idx,
      ref_fa = ref_fa,
      ref_fai = ref_fai,
      sample_name = sample_name,
      gq_binsize = gq_binsize,
      deepvariant_docker = deepvariant_docker
  }

  call call_variants {
    input:
      examples_tar = make_examples.examples_tar,
      sample_name = sample_name,
      deepvariant_model_path = deepvariant_model_path,
      deepvariant_docker = deepvariant_docker
  }

  call postprocess_variants {
    input:
      ref_fa = ref_fa,
      ref_fai = ref_fai,
      call_variants_output = call_variants.call_variants_output,
      gvcf_tfrecords_tar = make_examples.gvcf_tfrecords_tar,
      sample_name = sample_name,
      deepvariant_docker = deepvariant_docker
  }

  output {
      File vcf_gz = postprocess_variants.vcf_gz
      File gvcf_gz = postprocess_variants.gvcf_gz
  }
}


# DeepVariant make_examples
task make_examples {
  input {
    File bam
    File bam_idx
    File ref_fa
    File ref_fai
    String sample_name
    String deepvariant_docker

    Int? gq_binsize
  }

  Int disk_size = ceil( ( 2 * size([bam, ref_fa], "GB") ) + 20.0 )

  command <<<
    set -eu -o pipefail

    export SHELL=/bin/bash
    apt-get update -qq && apt-get install -y -qq samtools

    binsize_arg=""
    if [ -n "~{gq_binsize}" ]; then
        binsize_arg="--gvcf_gq_binsize ~{gq_binsize}"
    fi

    mkdir examples/ gvcf/ logs/
    output_fn="examples/~{sample_name}.tfrecord@$(nproc).gz"
    gvcf_fn="gvcf/~{sample_name}.gvcf.tfrecord@$(nproc).gz"

    seq 0 $(( `nproc` - 1 )) | NO_GCE_CHECK=True parallel --halt 2 -t --results logs/ \
      "/opt/deepvariant/bin/make_examples \
         --mode calling \
         --ref ~{ref_fa} \
         --reads ~{bam} \
         --examples '$output_fn' \
         --gvcf '$gvcf_fn' \
         --add_hp_channel=true \
         --task {} \
         $binsize_arg 2>&1"

    mkdir tar/
    tar -cvzf "tar/~{sample_name}.logs.tar.gz" -C logs/ .
    tar -cvf "tar/~{sample_name}.examples.tar" -C examples/ .
    tar -cvf "tar/~{sample_name}.gvcf_tfrecords.tar" -C gvcf/ .
  >>>

  runtime {
    docker: "~{deepvariant_docker}"
    cpu: "32"
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File examples_tar = glob("tar/*.examples.tar")[0]
    File gvcf_tfrecords_tar = glob("tar/*.gvcf_tfrecords.tar")[0]
    File logs_tar_gz = glob("tar/*.logs.tar.gz")[0]
  }
}


# DeepVariant call_variants
task call_variants {
  input {
    File examples_tar
    String sample_name
    String deepvariant_model_path
    String deepvariant_docker
  }

  Int disk_size = ceil( ( 2 * size(examples_tar, "GB") ) + 20.0 )

  parameter_meta {
    examples_tar: "stream"
  }

  command <<<
    set -eu -o pipefail

    mkdir examples output
    tar xvf "~{examples_tar}" -C examples/
    wait -n

    n_examples=$( find examples/ -type f -name "*.gz" | wc -l )

    NO_GCE_CHECK=True /opt/deepvariant/bin/call_variants \
      --outfile "output/~{sample_name}.call_variants.tfrecord.gz" \
      --examples "examples/~{sample_name}.tfrecord@$n_examples.gz" \
      --checkpoint ~{deepvariant_model_path}
  >>>

  runtime {
    docker: "~{deepvariant_docker}"
    cpu: "32"
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File call_variants_output = glob("output/*.gz")[0]
  }
}


# DeepVariant postprocess_variants
task postprocess_variants {
  input {
    File gvcf_tfrecords_tar
    File call_variants_output
    File ref_fa
    File ref_fai
    String sample_name
    String deepvariant_docker
  }

  Int disk_size = ceil( ( 2 * size([gvcf_tfrecords_tar, call_variants_output, ref_fa], "GB") ) + 20.0 )

  parameter_meta {
    ref_fa: "stream"
    gvcf_tfrecords_tar: "stream"
  }

  command <<<
    set -eu -o pipefail

    apt-get update -qq && apt-get install -y -qq samtools wget
    # download a multithreaded version of the tabix bgzip utility
    wget --quiet -O bgzip "https://github.com/dnanexus-rnd/GLnexus/blob/master/cli/dxapplet/resources/usr/local/bin/bgzip?raw=true"
    chmod +x bgzip

    mv ~{ref_fa} reference.fa
    mv ~{ref_fai} reference.fa.fai

    mkdir gvcf output
    tar xvf "~{gvcf_tfrecords_tar}" -C gvcf/
    n_gvcf_tfrecords=$(find gvcf/ -type f | wc -l)
    NO_GCE_CHECK=True /opt/deepvariant/bin/postprocess_variants \
      --ref reference.fa --infile "~{call_variants_output}" \
      --nonvariant_site_tfrecord_path "gvcf/~{sample_name}.gvcf.tfrecord@$n_gvcf_tfrecords.gz" \
      --outfile "output/~{sample_name}.vcf" \
      --gvcf_outfile "output/~{sample_name}.gvcf"

    ./bgzip -@ $(nproc) output/*.vcf &
    ./bgzip -@ $(nproc) output/*.gvcf
    wait -n
  >>>

  runtime {
    docker: "~{deepvariant_docker}"
    cpu: "8"
    preemptible: 3
    disks: "local-disk " + disk_size + " HDD"
  }

  output {
    File vcf_gz = glob("output/*.vcf.gz")[0]
    File gvcf_gz = glob("output/*.gvcf.gz")[0]
  }
}