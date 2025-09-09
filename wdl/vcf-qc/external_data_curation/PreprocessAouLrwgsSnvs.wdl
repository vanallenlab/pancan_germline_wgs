# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Preprocess existing AoU long-read WGS-based SNV/indel calls to match G2C QC expectations


version 1.0


import "QcTasks.wdl" as QcTasks


workflow PreprocessAouLrwgsSnvs {
  input {
	  File vcf = "gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_vcf/GRCh38/cohort_for_GLNexus_2023Q1_1027.g.vcf.bgz"
    File vcf_idx = "gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_vcf/GRCh38/cohort_for_GLNexus_2023Q1_1027.g.vcf.bgz.tbi"
    File? samples_list

    File ref_fasta
    File ref_fasta_idx
    Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
                             "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                             "chr13", "chr14", "chr15", "chr16", "chr17", 
                             "chr18", "chr19", "chr20", "chr21", "chr22", 
                             "chrX",  "chrY"]

    String g2c_pipeline_docker
  }

  call QcTasks.MakeHeaderFiller {}

  scatter ( contig in contigs ) {

    call CleanSnvs {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        contig = contig,
        ref_fasta = ref_fasta,
        ref_fasta_idx = ref_fasta_idx,
        supp_vcf_header = MakeHeaderFiller.supp_vcf_header,
        samples_list = samples_list,
        outfile_name = "AoU.lrWGS.snv_indel.cleaned." + contig + ".vcf.gz",
        g2c_pipeline_docker = g2c_pipeline_docker
    }

  }

  output {
    Array[File] cleaned_vcfs = CleanSnvs.cleaned_vcf
    Array[File] cleaned_vcf_idxs = CleanSnvs.cleaned_vcf_idx
  }
}


task CleanSnvs {
  input {
    File vcf
    File vcf_idx
    String contig

    File ref_fasta
    File ref_fasta_idx

    File supp_vcf_header

    File? samples_list

    String outfile_name
    
    Int? disk_gb
    Float mem_gb = 4
    Int n_cpu = 2

    String g2c_pipeline_docker    
  }

  String samples_cmd = if defined(samples_list) then "--samples-file ~{basename(select_first([samples_list]))} --force-samples" else ""

  Int default_disk_gb = ceil(2.5 * size(vcf, "GB")) + 10
  Int hdd_gb = select_first([disk_gb, default_disk_gb])
  Float sort_mem = 0.6 * mem_gb

  command <<<
    set -eu -o pipefail

    # Localize samples list to pwd if provided
    if [ ~{defined(samples_list)} ]; then
      cp ~{samples_list} ./
    fi

    # Filter & clean SNVs & indels
    bcftools view --regions "~{contig}" ~{samples_cmd} ~{vcf} \
    | bcftools norm -m -any -f ~{ref_fasta} -c s --threads ~{n_cpu} \
    | bcftools view -c 1 \
    | bcftools annotate -h ~{supp_vcf_header} -x "^FORMAT/GT" \
    | bcftools +fill-tags \
    | bcftools annotate --threads 2 \
      --set-id +'%CHROM\_%POS\_%REF\_%FIRST_ALT' \
      -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE,^FORMAT/GT" \
    | bcftools sort -m ~{sort_mem}G -Oz -o "~{outfile_name}"
    tabix -f -p vcf "~{outfile_name}"
  >>>

  output {
    File cleaned_vcf = "~{outfile_name}"
    File cleaned_vcf_idx = "~{outfile_name}.tbi"
  }

  runtime {
    docker: g2c_pipeline_docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + hdd_gb + " HDD"
    preemptible: 1
  }
}
