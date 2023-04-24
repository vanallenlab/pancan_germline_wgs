# The Genomic Architecture of Human Cancers
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Align PE WGS fastqs with bwa mem and preprocess with GATK4 best practices

# Originally based on a WDL by Tyler Chinsky <tyler_chinsky@dfci.harvard.edu>
# https://portal.firecloud.org/?return=terra#methods/tylerchinskydfci/paired_fastq_to_bam/5


version 1.0


# import the GATK pre-processing workflow
import "https://raw.githubusercontent.com/gatk-workflows/gatk4-data-processing/master/processing-for-variant-discovery-gatk4.wdl" as GATK


workflow WgsFastqToCram {
  input {
    String sample_name
    File fastq_1
    File fastq_2
    Boolean clean_fastqs = false
    Boolean cram_out = true

    String reference_name
    File reference_fasta
    File reference_fasta_index
    File reference_dict
    File? reference_alt
    File reference_sa 
    File reference_ann
    File reference_bwt
    File reference_pac
    File reference_amb
    File dbSNP_vcf
    File dbSNP_vcf_index
    File known_indels_sites_VCF
    File known_indels_sites_index

    String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
    String fastq_pair_docker = "vanallenlab/fastq-pair:latest"

    String fastq_suffix = ".fq.gz"
    Int clean_fastq_reads_per_hash_cell = 8
    Int clean_fastq_disk_scaling_factor = 6
    Float clean_fastq_mem_gb = 64
  }

  # Clean up FASTQs with fastq-pair
  if ( clean_fastqs ) {
    call CleanFastqs {
      input:
        fastq_1 = fastq_1,
        fastq_2 = fastq_2,
        fastq_suffix = fastq_suffix,
        reads_per_hash_cell = clean_fastq_reads_per_hash_cell,
        disk_scaling_factor = clean_fastq_disk_scaling_factor,
        mem_gb = clean_fastq_mem_gb,
        docker = fastq_pair_docker
    }
  }
  File fq1 = select_first([CleanFastqs.fq_paired_1, fastq_1])
  File fq2 = select_first([CleanFastqs.fq_paired_2, fastq_2])
  File fq_singles = select_first([CleanFastqs.fq_singles, []])

  # Align with BWA-MEM
  call Bwa as AlignPairs {}
  scatter ( fastq in fq_singles ) {
    call Bwa as AlignSingles {}
  }

  # # Convert the uBAM to bam using the GATK pre-processing workflow
  # call GATK.PreProcessingForVariantDiscovery_GATK4 {
  #   input:
  #     sample_name = sample_name,
  #     ref_name = reference_name,

  #     flowcell_unmapped_bams_list = WriteUBAMPath.ubam_path_file,
  #     unmapped_bam_suffix = ".unmapped.bam",
    
  #     ref_fasta = reference_fasta,
  #     ref_fasta_index = reference_fasta_index,
  #     ref_dict = reference_dict,
  #     ref_alt = reference_alt,
  #     ref_sa  = reference_sa,
  #     ref_ann = reference_ann,
  #     ref_bwt = reference_bwt,
  #     ref_pac = reference_pac,
  #     ref_amb = reference_amb,
  #     dbSNP_vcf = dbSNP_vcf,
  #     dbSNP_vcf_index = dbSNP_vcf_index,
  #     known_indels_sites_VCFs = [known_indels_sites_VCF],
  #     known_indels_sites_indices = [known_indels_sites_index],

  #     bwa_commandline = "bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta",
  #     compression_level = 5,
    
  #     gatk_docker = "us.gcr.io/broad-gatk/gatk:4.2.0.0",
  #     gatk_path = "/gatk/gatk",
  #     gotc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710",
  #     gotc_path = "/usr/gitc/",
  #     python_docker = "python:2.7",

  #     flowcell_small_disk = 500,
  #     flowcell_medium_disk = 1000,
  #     agg_small_disk = 500,
  #     agg_medium_disk = 1000,
  #     agg_large_disk = 1500,

  #     preemptible_tries = 3
  # }

  # if ( cram_out ) {
  #   call Bam2Cram {
  #     input:
  #       bam = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam,
  #       bai = PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index,
  #       ref_fa = reference_fasta,
  #       ref_fai = reference_fasta_index,
  #       docker = gatk_docker
  #   }
  # }

  # output {
  #   File cram_or_bam = select_first([Bam2Cram.cram, PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam])
  #   File cram_or_bam_index = select_first([Bam2Cram.crai, PreProcessingForVariantDiscovery_GATK4.analysis_ready_bam_index])
  # }
}

# Clean up FASTQs with fastq-pair
task CleanFastqs {
  input {
    File fastq_1
    File fastq_2
    String docker
    String fastq_suffix = ".fq.gz"
    Int reads_per_hash_cell = 3
    Int disk_scaling_factor = 10
    Float mem_gb = 31.5
  }

  String f1_basename = basename(fastq_1, fastq_suffix)
  String f2_basename = basename(fastq_2, fastq_suffix)

  Int disk_gb_base = ceil(3 * size([fastq_1, fastq_2], "GB"))
  Int disk_gb = (disk_scaling_factor * disk_gb_base) + 10

  command <<<
    set -eu -o pipefail

    # log resource usage for debugging purposes
    function runtimeInfo() {
      echo [$(date)]
      echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
      echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
      echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
    }
    while true; do
      runtimeInfo
      ls -lh
      sleep 60
    done &

    # Check whether fastqs are gzipped and uncompress if necessary
    if [ $( file ~{fastq_1} | fgrep gzip | wc -l ) -gt 0 ]; then
      zcat ~{fastq_1} > ~{f1_basename}.fastq
    else
      mv ~{fastq_1} > ~{f1_basename}.fastq
    fi
    if [ $( file ~{fastq_2} | fgrep gzip | wc -l ) -gt 0 ]; then
      zcat ~{fastq_1} > ~{f2_basename}.fastq
    else
      mv ~{fastq_1} > ~{f2_basename}.fastq
    fi

    # Get number of reads
    n_reads=$( cat ~{f1_basename}.fastq | wc -l | awk '{ print $1 / 4 }' | cut -f1 -d\. )
    echo -e "\nCounted $n_reads total reads\n"

    # Set hash size to n_reads / reads_per_hash_cell
    hash_size=$( echo $n_reads | awk -v denom=~{reads_per_hash_cell} '{ printf "%.0f\n", $1 / denom }' | cut -f1 -d\. )
    echo -e "\nSetting hash size to $hash_size\n"

    fastq_pair -t $hash_size ~{f1_basename}.fastq ~{f2_basename}.fastq
    for out in paired single; do
      gzip -f ~{f1_basename}.fastq.$out.fq
      gzip -f ~{f2_basename}.fastq.$out.fq
    done
  >>>

  output {
    File fq_paired_1 = "~{f1_basename}.fastq.paired.fq.gz"
    File fq_paired_2 = "~{f2_basename}.fastq.paired.fq.gz"
    Array[File] fq_singles = ["~{f1_basename}.fastq.single.fq.gz",
                              "~{f2_basename}.fastq.single.fq.gz"]
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}

# Write paired fastqs to uBAM
task PairedFastqsToUnmappedBAM {
  input {
    String sample_name
    File fastq_1
    File fastq_2
    String readgroup_name

    String gatk_docker
    String gatk_path

    # Runtime parameters
    Int addtional_disk_space_gb = 20
    Int preemptible_attempts = 3
    Int machine_mem_gb = 4
  }

  Int command_mem_gb = machine_mem_gb - 1
  Int disk_space_gb = ceil(6 * size([fastq_1, fastq_2], "GB")) + addtional_disk_space_gb

  command <<<
    set -euo pipefail

    # log resource usage for debugging purposes
    function runtimeInfo() {
      echo [$(date)]
      echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
      echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
      echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
    }
    while true; do
      runtimeInfo
      sleep 60s
    done &

    ~{gatk_path} \
    --java-options "-Xmx~{command_mem_gb}g" \
    FastqToSam \
    --FASTQ ~{fastq_1} \
    --FASTQ2 ~{fastq_2} \
    --OUTPUT ~{sample_name}.unmapped.bam \
    --COMPRESSION_LEVEL 5 \
    --READ_GROUP_NAME ~{readgroup_name} \
    --SAMPLE_NAME ~{sample_name}
  >>>

  output {
    File ubam = "~{sample_name}.unmapped.bam"
  }

  runtime {
    docker: gatk_docker
    memory: machine_mem_gb + " GB"
    cpu: 2
    disks: "local-disk " + disk_space_gb + " HDD"
    preemptible: preemptible_attempts
  }
}

# Write uBAM gs:// path to file (required by GATK preprocessing)
task WriteUBAMPath {
  input {
    String ubam
    String outfile
  }

  command <<<
    echo ~{ubam} > "~{outfile}"
  >>>

  output {
    File ubam_path_file = "~{outfile}"
  }
  runtime {
    cpu: 1
    memory: "2 GB"
    docker: "ubuntu:latest"
    preemptible: 3
  }
}

# Convert BAM to CRAM
task Bam2Cram {
  input {
    File bam
    File bai
    File ref_fa
    File ref_fai

    String docker
  }

  String cram_filename = basename(bam, ".bam") + ".cram"

  Int disk_gb = ceil(3 * size(bam, "GB")) + 10

  command <<<
    set -euo pipefail

    # log resource usage for debugging purposes
    function runtimeInfo() {
      echo [$(date)]
      echo \* CPU usage: $(top -bn 2 -d 0.01 | grep '^%Cpu' | tail -n 1 | awk '{print $2}')%
      echo \* Memory usage: $(free -m | grep Mem | awk '{ OFMT="%.0f"; print ($3/$2)*100; }')%
      echo \* Disk usage: $(df | grep cromwell_root | awk '{ print $5 }')
    }
    while true; do
      runtimeInfo
      sleep 60s
    done &

    samtools view \
      -C \
      -T ~{ref_fa} \
      -o ~{cram_filename} \
      --threads 4 \
      ~{bam}

    samtools index \
      -@ 4 \
      ~{cram_filename}
  >>>

  output {
    File cram = "~{cram_filename}"
    File crai = "~{cram_filename}.crai"
  }

  runtime {
    docker: docker
    memory: "7.5 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}
