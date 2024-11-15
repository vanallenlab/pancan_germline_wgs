# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Remove a list of predefined outlier samples from the inputs to GATK-SV module 05 (ClusterBatch.wdl)


version 1.0


workflow ExcludeClusteredOutliers {
  input {
    File del_bed
    File dup_bed
    File manta_vcf_tar
    File melt_vcf_tar
    File wham_vcf_tar
    File exclude_samples_list
    String linux_docker
  }

  # Exclude samples from depth-based CNV BED files
  call ExcludeOutliersFromDepthBed as CleanDelBed {
    input:
      bed = del_bed,
      exclude_samples_list = exclude_samples_list,
      docker = linux_docker
  }
  call ExcludeOutliersFromDepthBed as CleanDupBed {
    input:
      bed = dup_bed,
      exclude_samples_list = exclude_samples_list,
      docker = linux_docker
  }

  # Exclude samples from PE/SR-based SV VCF tarballs
  call ExcludeOutliersFromVcfTarball as CleanManta {
    input:
      tarball = manta_vcf_tar,
      exclude_samples_list = exclude_samples_list,
      docker = linux_docker
  }
  call ExcludeOutliersFromVcfTarball as CleanMelt {
    input:
      tarball = melt_vcf_tar,
      exclude_samples_list = exclude_samples_list,
      docker = linux_docker
  }
  call ExcludeOutliersFromVcfTarball as CleanWham {
    input:
      tarball = wham_vcf_tar,
      exclude_samples_list = exclude_samples_list,
      docker = linux_docker
  }

  output {
    File del_bed_cleaned = CleanDelBed.cleaned_bed
    File dup_bed_cleaned = CleanDupBed.cleaned_bed
    File manta_vcf_tar_cleaned = CleanManta.cleaned_tarball
    File melt_vcf_tar_cleaned = CleanMelt.cleaned_tarball
    File wham_vcf_tar_cleaned = CleanWham.cleaned_tarball
  }
}


task ExcludeOutliersFromDepthBed {
  input {
    File bed
    File exclude_samples_list
    String docker
  }

  String output_bed = basename(bed, ".bed.gz") + ".no_outliers.bed.gz"

  Int disk_gb = ceil(3 * size(bed, "GB")) + 10

  command <<<
    set -eu -o pipefail

    zcat ~{bed} | fgrep -wvf ~{exclude_samples_list} | bgzip -c > ~{output_bed}
  >>>

  output {
    File cleaned_bed = "~{output_bed}"
  }

  runtime {
    docker: docker
    memory: "2.0 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task ExcludeOutliersFromVcfTarball {
  input {
    File tarball
    File exclude_samples_list
    String docker
  }

  String output_tarball = basename(tarball, ".tar.gz") + ".no_outliers.tar.gz"

  Int disk_gb = ceil(10 * size(tarball, "GB")) + 10

  command <<<
    set -eu -o pipefail

    mkdir tar_contents && cd tar_contents

    tar -xzvf ~{tarball}

    while read sid; do
      rm *.$sid.vcf.gz
    done < ~{exclude_samples_list}

    tar -czvf ~{output_tarball} ./*
  >>>

  output {
    File cleaned_tarball = "~{output_tarball}"
  }

  runtime {
    docker: docker
    memory: "2.0 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}
