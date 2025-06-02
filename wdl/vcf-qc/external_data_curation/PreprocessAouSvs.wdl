# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Preprocess existing AoU SV VCFs for short- and long-read WGS to match G2C QC expectations


version 1.0


import "QcTasks.wdl" as QcTasks


workflow PreprocessAouSvs {
  input {
    String srwgs_vcf_prefix = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/structural_variants/vcf/full/AoU_srWGS_SV.v8"
    File lrwgs_vcf = "gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_sv/GRCh38/integrated_sv_with_hprc_year_1_more_stringent.vcf.gz"
    File lrwgs_vcf_idx = "gs://fc-aou-datasets-controlled/v7/wgs/long_read/joint_sv/GRCh38/integrated_sv_with_hprc_year_1_more_stringent.vcf.gz.tbi"

    File? samples_list

    Array[String] contigs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
                             "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
                             "chr13", "chr14", "chr15", "chr16", "chr17", 
                             "chr18", "chr19", "chr20", "chr21", "chr22", 
                             "chrX",  "chrY"]

    String g2c_pipeline_docker
  }

  call QcTasks.MakeHeaderFiller {}

  scatter ( contig in contigs ) {
    # Curate srWGS SVs
    call CurateSrwgsSvs {
      input:
        vcf = srwgs_vcf_prefix + "." + contig + ".vcf.gz",
        vcf_idx = srwgs_vcf_prefix + "." + contig + ".vcf.gz.tbi",
        out_vcf_fname = "AoU.srWGS.sv.cleaned." + contig + ".vcf.gz",
        samples_list = samples_list,
        supp_vcf_header = MakeHeaderFiller.supp_vcf_header,
        g2c_pipeline_docker = g2c_pipeline_docker
    }

    # Curate lrWGS SVs
    call CurateLrwgsSvs {
      input:
        vcf = lrwgs_vcf,
        vcf_idx = lrwgs_vcf_idx,
        out_vcf_fname = "AoU.lrWGS.sv.cleaned." + contig + ".vcf.gz",
        contig = contig,
        samples_list = samples_list,
        supp_vcf_header = MakeHeaderFiller.supp_vcf_header,
        g2c_pipeline_docker = g2c_pipeline_docker
    }
  }

  output {
    Array[File] cleaned_srwgs_vcfs = CurateSrwgsSvs.vcf_out
    Array[File] cleaned_srwgs_vcf_idxs = CurateSrwgsSvs.idx_out
    Array[File] cleaned_lrwgs_vcfs = CurateLrwgsSvs.vcf_out
    Array[File] cleaned_lrwgs_vcf_idxs = CurateLrwgsSvs.idx_out
  }
}


task CurateSrwgsSvs {
  input {
    File vcf
    File vcf_idx
    String out_vcf_fname
    File? samples_list
    File supp_vcf_header
    String g2c_pipeline_docker
  }

  String out_tbi_fname = out_vcf_fname + ".tbi"
  Int disk_gb = ceil(3 * size(vcf, "GB")) + 10
  String sample_subset_cmd = if defined(samples_list) then "--samples-file ~{basename(select_first([samples_list]))} --force-samples" else ""

  command <<<
    set -eu -o pipefail

    # Localize samples list to pwd if provided
    if [ ~{defined(samples_list)} ]; then
      cp ~{samples_list} ./
    fi

    # Filter & clean VCF
    bcftools view -f PASS,MULTIALLELIC ~{sample_subset_cmd} ~{vcf} \
    | bcftools +fill-tags \
    | bcftools annotate -h ~{supp_vcf_header} \
    | bcftools annotate --threads 2 --include 'AC>0 & AN>0' \
      -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE,^FORMAT/GT,FORMAT/RD_CN,^FILTER/PASS,FILTER/MULTIALLELIC" \
      -Oz -o "~{out_vcf_fname}"
    tabix -f -p vcf "~{out_vcf_fname}"
  >>>

  output {
    File vcf_out = "~{out_vcf_fname}"
    File idx_out = "~{out_tbi_fname}"
  }

  runtime {
    docker: g2c_pipeline_docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}


task CurateLrwgsSvs {
  input {
    File vcf
    File vcf_idx
    String out_vcf_fname
    String contig
    File? samples_list
    File supp_vcf_header
    String g2c_pipeline_docker
  }

  String out_tbi_fname = out_vcf_fname + ".tbi"
  Int disk_gb = ceil(3 * size(vcf, "GB")) + 10
  String sample_subset_cmd = if defined(samples_list) then "--samples-file ~{basename(select_first([samples_list]))} --force-samples" else ""

  command <<<
    set -eu -o pipefail

    # Localize samples list to pwd if provided
    if [ ~{defined(samples_list)} ]; then
      cp ~{samples_list} ./
    fi

    # Filter & clean VCF
    bcftools view --regions "~{contig}" ~{sample_subset_cmd} ~{vcf} \
    | bcftools annotate -h ~{supp_vcf_header} -x "^INFO/SVTYPE,INFO/SVLEN,^FORMAT/GT" \
    | /opt/pancan_germline_wgs/scripts/data_management/external_data_curation/postprocess_hgsvc_vcfs.py stdin stdout --sv --no-gt \
    | bcftools +fill-tags \
    | bcftools view -c 1 \
    | bcftools annotate --threads 2 \
      -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AN,INFO/AC,INFO/AF,INFO/CN_NONREF_COUNT,INFO/CN_NONREF_FREQ,INFO/AC_Het,INFO/AC_Hom,INFO/AC_Hemi,INFO/HWE,^FORMAT/GT" \
      -Oz -o "~{out_vcf_fname}"
    tabix -f -p vcf "~{out_vcf_fname}"
  >>>

  output {
    File vcf_out = "~{out_vcf_fname}"
    File idx_out = "~{out_tbi_fname}"
  }

  runtime {
    docker: g2c_pipeline_docker
    memory: "4 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}
