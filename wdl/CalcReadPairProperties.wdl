# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Compute read length and insert size for a single short-read WGS library

# Uses "GetMultipleMetrics" from GATK-SV MELT.wdl, which is reproduced in this
# WDL rather than using an import statement to ensure cross-platform portability
# (e.g., Terra import vs. All of Us RW dependencies.zip)


version 1.0


workflow CalcReadPairProperties {
  input {
    String sample_id
    File bam_or_cram_file
    File bam_or_cram_index
    File reference_fasta
    File reference_index

    Boolean stage_output = false
    String? cohort
    String staging_bucket_base = "gs://dfci-g2c-inputs/"

    String gatk_docker
    String gcloud_docker = "google/cloud-sdk:latest"
  }

  call GetMultipleMetrics {
    input:
      sample_id = sample_id,
      bam_or_cram_file = bam_or_cram_file,
      bam_or_cram_index = bam_or_cram_index,
      reference_fasta = reference_fasta,
      reference_index = reference_index,
      gatk_docker = gatk_docker
  }

  if ( stage_output ) {
    call StageOutput {
      input:
        files_to_stage = [GetMultipleMetrics.simple_read_metrics],
        cohort = cohort,
        staging_bucket_base = staging_bucket_base,
        gcloud_docker = gcloud_docker
    }
  }

  output {
    Array[File] metrics_files = GetMultipleMetrics.metrics_files
    File metrics_table_file = GetMultipleMetrics.metrics_table_file
  }
}


task GetMultipleMetrics {
  input {
    File bam_or_cram_file
    File bam_or_cram_index
    File reference_fasta
    File reference_index
    String sample_id
    String gatk_docker
  }

  String metrics_base = "multiple_metrics"
  String alignment_summary_filename = metrics_base + ".alignment_summary_metrics"
  String insert_size_filename = metrics_base + ".insert_size_metrics"
  String gc_bias_filename = metrics_base + ".gc_bias.summary_metrics"
  String metrics_table_filename = metrics_base + "_table.tsv"
  String g2c_metrics_filename = sample_id + ".read_metrics.tsv"

  Float bam_or_cram_size = size(bam_or_cram_file, "GiB")
  Float disk_overhead = 20.0
  Float ref_size = size(reference_fasta, "GiB")
  Int vm_disk_size = ceil(bam_or_cram_size + ref_size + disk_overhead)

  Int cpu_cores = 2
  Float mem_use_gb = 6.0
  Float java_mem_pad_gb = cpu_cores * 1.0
  Float mem_size_gb = mem_use_gb + java_mem_pad_gb

  command <<<

    set -Eeuo pipefail

    gatk --java-options -Xmx~{java_mem_mb}m CollectMultipleMetrics \
      -I "~{bam_or_cram_file}" \
      -O "~{metrics_base}" \
      -R "~{reference_fasta}" \
      --ASSUME_SORTED true \
      --PROGRAM null \
      --PROGRAM CollectAlignmentSummaryMetrics \
      --PROGRAM CollectInsertSizeMetrics \
      --METRIC_ACCUMULATION_LEVEL null \
      --METRIC_ACCUMULATION_LEVEL SAMPLE

    function transpose_table() {
          cat \
          | awk ' {
              for (col = 1; col <= NF; ++col) {
                table[NR, col] = $col
              }
              if(NF > num_cols) {
                num_cols = NF
              }
            } END {
              for (row = 1; row <= num_cols; ++row) {
                printf "%s", table[1, row]
                for (col = 2; col <= NR; ++col) {
                  printf "\t%s", table[col, row]
                }
                printf "\n"
              }
            }'
        }

    # get alignment summary metrics for all PAIR from sample and transpose to tsv
    grep -A4 "^## METRICS CLASS" "~{alignment_summary_filename}" \
      | grep -E "^CATEGORY\s|^PAIR\s" \
      | transpose_table \
      | grep -Ev "^CATEGORY\s|^LIBRARY\s|^SAMPLE\s|^READ_GROUP\s" \
      > "~{metrics_table_filename}"

    # get insert size metrics from sample and transpose to tsv
    grep -A1 MEDIAN_INSERT_SIZE "~{insert_size_filename}"  \
      | transpose_table \
      | grep -Ev "^LIBRARY\s|^SAMPLE\s|^READ_GROUP\s" \
      >> "~{metrics_table_filename}"

    # Format values for G2C
    rl=$( fgrep -w "MEAN_READ_LENGTH" "~{metrics_table_filename}" | cut -f2 )
    isize=$( fgrep -w "MEAN_INSERT_SIZE" "~{metrics_table_filename}" | cut -f2 )
    echo -e "#Sample\tinsert_size\tread_length" > "~{g2c_metrics_filename}"
    echo -e "~{sample_id}\t$isize\t$rl" >> "~{g2c_metrics_filename}"
  >>>

  runtime {
    cpu: cpu_cores
    memory: mem_use_gb + " GiB"
    disks: "local-disk " + vm_disk_size + " HDD"
    bootDiskSizeGb: 10
    docker: gatk_docker
    preemptible: 1
    maxRetries: 1
  }
  Int java_mem_mb = round(1024 * (7.0 - java_mem_pad_gb))

  output {
    Array[File] metrics_files = glob("~{metrics_base}.*")
    File metrics_table_file = metrics_table_filename
    File simple_read_metrics = g2c_metrics_filename
  }
}


task StageOutput {
  input {
    Array[File] files_to_stage
    String? cohort
    String staging_bucket_base
    String gcloud_docker
  }

  String target_bucket = staging_bucket_base + select_first([cohort]) + "/gatk-sv/metrics/"

  Int disk_gb = ceil(size(files_to_stage, "GB") + 20)

  command <<<
    set -eu -o pipefail

    gsutil -m cp ~{sep=" " files_to_stage} ~{target_bucket}
  >>>

  runtime {
    cpu: 1
    memory: "2.0 GiB"
    disks: "local-disk " + disk_gb + " HDD"
    bootDiskSizeGb: 10
    docker: gcloud_docker
    preemptible: 3
    maxRetries: 1
  }

  output {}
}
