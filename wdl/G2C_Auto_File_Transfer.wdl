# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

#The purpose of this code is to transfer files from Terra to Google Cloud Storage.


version 1.0


task gatk_hc_reblocked{
  input {
    #Go into /cohort/gatk-hc/reblocked
    File reblocked_gvcf
    File reblocked_gvcf_idx
    String sample_id
    String cohort
  }
  Int disk_gb = ceil(1.5 * size(reblocked_gvcf, "GB")) + 10
  String reblocked_gvcf_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-hc/reblocked/~{sample_id}.reblocked.g.vcf.gz"
  String reblocked_gvcf_idx_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-hc/reblocked/~{sample_id}.reblocked.g.vcf.gz.tbi"
  command <<<
    gsutil -m cp ~{reblocked_gvcf} ~{reblocked_gvcf_pathway}
    gsutil -m cp ~{reblocked_gvcf_idx} ~{reblocked_gvcf_idx_pathway}
  >>>
  runtime{
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    disks: "local-disk " + disk_gb + " HDD"
  }

}

task gatk_sv_coverage{
  input {
    #Go into /cohort/gatk-sv/coverage
    File coverage_counts
    String sample_id
    String cohort
  }
  String coverage_counts_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/coverage/~{sample_id}.counts.tsv.gz"
  
  command <<<
    gsutil -m cp ~{coverage_counts} ~{coverage_counts_pathway}
  >>>
  runtime{
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
  }
}

task gatk_sv_metrics{
  input {
    #Go into /cohort/gatk-sv/metrics
    Array[File] sample_metric_files
    String sample_id
    String cohort
  }
  String pe_file_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/~{sample_id}.pe-file.tsv"
  String raw_counts_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/~{sample_id}.raw-counts.tsv"
  String sr_file_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/~{sample_id}.sr-file.tsv"
  String manta_vcf_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/manta_~{sample_id}.vcf.tsv"
  String melt_vcf_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/melt_~{sample_id}.vcf.tsv"
  String wham_vcf_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/wham_~{sample_id}.vcf.tsv"
  
  command <<<
    gsutil -m cp ~{sep=" " sample_metric_files} "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/"

  >>>
  runtime{
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
  }
}

task gatk_sv_pesr{
  input {
    #Go into /cohort/gatk-sv/pesr/
    File pesr_disc
    File pesr_disc_index
    File pesr_sd
    File pesr_sd_index
    File pesr_split
    File pesr_split_index
    String sample_id
    String cohort
  }
  String pesr_disc_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/pesr/~{sample_id}.pe.txt.gz"
  String pesr_disc_index_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/pesr/~{sample_id}.pe.txt.gz.tbi"
  
  String pesr_sd_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/pesr/~{sample_id}.sd.txt.gz"
  String pesr_sd_index_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/pesr/~{sample_id}.sd.txt.gz.tbi"
  
  String pesr_split_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/pesr/~{sample_id}.sr.txt.gz"
  String pesr_split_index_pathway = "gs://dfci-g2c-inputs/~{cohort}/gatk-sv/pesr/~{sample_id}.sr.txt.gz.tbi"
  
  command <<<
    gsutil -m cp ~{pesr_disc} ~{pesr_disc_pathway}
    gsutil -m cp ~{pesr_disc_index} ~{pesr_disc_index_pathway}
    
    gsutil -m cp ~{pesr_sd} ~{pesr_sd_pathway}
    gsutil -m cp ~{pesr_sd_index} ~{pesr_sd_index_pathway}
    
    gsutil -m cp ~{pesr_split} ~{pesr_split_pathway}
    gsutil -m cp ~{pesr_split_index} ~{pesr_split_index_pathway}
  >>>
  runtime{
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
  }

}
task manta{
  input {
    #Go into /cohort/manta
    File manta_vcf
    File manta_vcf_tbi
    String sample_id
    String cohort
  }
  Int disk_gb = ceil(1.5 * size(manta_vcf, "GB")) + 10
  String manta_vcf_pathway = "gs://dfci-g2c-inputs/~{cohort}/manta/~{sample_id}.manta.vcf.gz"
  String manta_vcf_tbi_pathway = "gs://dfci-g2c-inputs/~{cohort}/manta/~{sample_id}.manta.vcf.gz.tbi"
  
  command <<<
    gsutil -m cp ~{manta_vcf} ~{manta_vcf_pathway}
    gsutil -m cp ~{manta_vcf_tbi} ~{manta_vcf_tbi_pathway}
  >>>
  runtime{
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    disks: "local-disk " + disk_gb + " HDD"
  }

}

task melt{
  input {
    #Go into /cohort/melt
    File melt_vcf
    File melt_vcf_tbi
    String sample_id
    String cohort

  }
  Int disk_gb = ceil(1.5 * size(melt_vcf, "GB")) + 10
  String melt_vcf_pathway = "gs://dfci-g2c-inputs/~{cohort}/melt/~{sample_id}.melt.vcf.gz"
  String melt_vcf_tbi_pathway = "gs://dfci-g2c-inputs/~{cohort}/melt/~{sample_id}.melt.vcf.gz.tbi"
  
  command <<<
    gsutil -m cp ~{melt_vcf} ~{melt_vcf_pathway}
    gsutil -m cp ~{melt_vcf_tbi} ~{melt_vcf_tbi_pathway}
  >>>
  runtime{
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    disks: "local-disk " + disk_gb + " HDD"
  }
}

task wham{
  input {
    #Go into /cohort/wham
    File wham_vcf
    File wham_index
    String sample_id
    String cohort
  }
  Int disk_gb = ceil(1.5 * size(wham_vcf, "GB")) + 10
  #Should these be Strings or File variables
  String wham_vcf_pathway = "gs://dfci-g2c-inputs/~{cohort}/wham/~{sample_id}.wham.vcf.gz"
  String wham_vcf_index_pathway = "gs://dfci-g2c-inputs/~{cohort}/wham/~{sample_id}.wham.vcf.gz.tbi"
  
  command <<<
    gsutil -m cp ~{wham_vcf} ~{wham_vcf_pathway}
    gsutil -m cp ~{wham_index} ~{wham_vcf_index_pathway}
  >>>
  
  runtime{
    docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    disks: "local-disk " + disk_gb + " HDD"
  }
}


# Define workflow
workflow store_in_gcs {

  input {
  	String cohort
  	String sample_id
  }
  
  call gatk_hc{
    input:
      sample_id = sample_id,
      cohort = cohort
  }
  call gatk_hc_reblocked{
    input:
      sample_id = sample_id,
      cohort = cohort
  }
  call gatk_sv_coverage{
  	input:
      sample_id = sample_id,
      cohort = cohort
  }
  call gatk_sv_metrics{
  	input:
      sample_id = sample_id,
      cohort = cohort
  }
  call gatk_sv_pesr{
  	input:
      sample_id = sample_id,
      cohort = cohort
  }
  call manta{
  	input:
      sample_id = sample_id,
      cohort = cohort
  }
  call melt{
 	input:
      sample_id = sample_id,
      cohort = cohort
  }
  call wham{
  	input:
      sample_id = sample_id,
      cohort = cohort
  }
}