# The Genomic Architecture of Human Cancers
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Align PE WGS fastqs with bwa mem and preprocess with GATK4 best practices

# Originally based on a WDL by Tyler Chinsky <tyler_chinsky@dfci.harvard.edu>
# https://portal.firecloud.org/?return=terra#methods/tylerchinskydfci/paired_fastq_to_bam/5


version 1.0


# import GATK pre-processing workflow
import "https://raw.githubusercontent.com/gatk-workflows/gatk4-data-processing/master/processing-for-variant-discovery-gatk4.wdl" as GATK


workflow WgsFastqToCram {
  input {
    String sample_name
    File fastq_1
    File fastq_2
    Boolean clean_fastqs = false
    Int n_fastq_shards = 1
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

    String fastq_pair_docker = "vanallenlab/fastq-pair:latest"
    String fastqsplitter_docker = "vanallenlab/fastqsplitter:latest"
    String gotc_docker = "us.gcr.io/broad-gotc-prod/genomes-in-the-cloud:2.4.7-1603303710"
    String gatk_docker = "broadinstitute/gatk:latest"
    String gatk_path = "/gatk/gatk"
    String python_docker = "python:2.7"

    String fastq_suffix = ".fq.gz"
    Int clean_fastq_reads_per_hash_cell = 8
    Int clean_fastq_disk_scaling_factor = 6

    Float clean_fastq_mem_gb = 64
    Int? sort_and_fix_tags_disk_size
    Float? bwa_mem_gb
    Int? bwa_num_cpu
    Int? bwa_disk_size
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
  Array[File] fq_singles = select_first([CleanFastqs.fq_singles, []])

  # Split paired fastqs to parallelize alignment
  if ( n_fastq_shards > 1 ) {
    call SplitFastq as SplitFq1 {
      input:
        fastq = fq1,
        n_shards = n_fastq_shards,
        docker = fastqsplitter_docker
    }
    call SplitFastq as SplitFq2 {
      input:
        fastq = fq2,
        n_shards = n_fastq_shards,
        docker = fastqsplitter_docker
    }
  }
  Array[File] sharded_fq1 = select_first([SplitFq1.fq_shards, [fq1]])
  Array[File] sharded_fq2 = select_first([SplitFq2.fq_shards, [fq2]])
  Array[Pair[File, File]] sharded_fq_pairs = zip(sharded_fq1, sharded_fq2)

  # Align paired fastqs with BWA
  scatter ( fq_pair in sharded_fq_pairs ) {
    call Bwa as AlignPairs {
      input:
        input_fastqs = [fq_pair.left, fq_pair.right],
        output_basename = sample_name + "." + basename(fq_pair.left, ".fq.gz"),
        ref_fasta = reference_fasta,
        ref_fasta_index = reference_fasta_index,
        ref_dict = reference_dict,
        ref_alt = reference_alt,
        ref_amb = reference_amb,
        ref_ann = reference_ann,
        ref_bwt = reference_bwt,
        ref_pac = reference_pac,
        ref_sa = reference_sa,
        mem_gb = bwa_mem_gb,
        num_cpu = bwa_num_cpu,
        disk_size = bwa_disk_size,
        docker = gotc_docker
    }
  }

  scatter ( fastq in fq_singles ) {
    call Bwa as AlignSingles {
      input:
        input_fastqs = [fastq],
        output_basename = sample_name,
        ref_fasta = reference_fasta,
        ref_fasta_index = reference_fasta_index,
        ref_dict = reference_dict,
        ref_alt = reference_alt,
        ref_amb = reference_amb,
        ref_ann = reference_ann,
        ref_bwt = reference_bwt,
        ref_pac = reference_pac,
        ref_sa = reference_sa,
        mem_gb = bwa_mem_gb,
        num_cpu = bwa_num_cpu,
        disk_size = bwa_disk_size,
        docker = gotc_docker
    }
  }

  # Merge BAMs, mark duplicates
  Array[File] bams_for_gatk = flatten([AlignPairs.output_bam, AlignSingles.output_bam])
  call GATK.MarkDuplicates {
    input:
      input_bams = bams_for_gatk,
      output_bam_basename = sample_name + ".aligned.unsorted.duplicates_marked",
      metrics_filename = sample_name + ".duplicate_metrics",
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      disk_size = ceil(3 * size(bams_for_gatk, "GB")) + 10,
      compression_level = 2,
      preemptible_tries = 5
  }

  # Sort deduped BAM file and fix tags
  call SortAndFixTags {
    input:
      input_bam = MarkDuplicates.output_bam,
      output_bam_basename = sample_name + ".aligned.duplicate_marked.sorted",
      ref_dict = reference_dict,
      ref_fasta = reference_fasta,
      ref_fasta_index = reference_fasta_index,
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      disk_size = select_first([sort_and_fix_tags_disk_size, ceil(4 * size(MarkDuplicates.output_bam, "GB")) + 10]),
      preemptible_tries = 5,
      compression_level = 2
  }

  # Add RG tag (required for BaseRecalibrator)
  call AddRg {
    input:
      bam = SortAndFixTags.output_bam,
      bam_index = SortAndFixTags.output_bam_index,
      sample_name = sample_name,
      gatk_path = gatk_path,
      docker = gatk_docker
  }

  # Create list of sequences for scatter-gather parallelization 
  call GATK.CreateSequenceGroupingTSV {
    input:
      ref_dict = reference_dict,
      docker_image = python_docker,
      preemptible_tries = 5
  }
  
  # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM in parallel
  scatter ( subgroup in CreateSequenceGroupingTSV.sequence_grouping ) {
    # Generate the recalibration model by interval
    call GATK.BaseRecalibrator {
      input:
        input_bam = AddRg.bam_wRG,
        input_bam_index = AddRg.bam_wRG_idx,
        recalibration_report_filename = sample_name + ".recal_data.csv",
        sequence_group_interval = subgroup,
        dbSNP_vcf = dbSNP_vcf,
        dbSNP_vcf_index = dbSNP_vcf_index,
        known_indels_sites_VCFs = [known_indels_sites_VCF],
        known_indels_sites_indices = [known_indels_sites_index],
        ref_dict = reference_dict,
        ref_fasta = reference_fasta,
        ref_fasta_index = reference_fasta_index,
        docker_image = gatk_docker,
        gatk_path = gatk_path,
        disk_size = ceil(2 * size(AddRg.bam_wRG, "GB")) + 10,
        preemptible_tries = 5
    }  
  }  
  
  # Merge the recalibration reports resulting from by-interval recalibration
  call GATK.GatherBqsrReports {
    input:
      input_bqsr_reports = BaseRecalibrator.recalibration_report,
      output_report_filename = sample_name + ".recal_data.csv",
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      disk_size = 100,
      preemptible_tries = 5
  }

  scatter ( subgroup in CreateSequenceGroupingTSV.sequence_grouping_with_unmapped ) {

    # Apply the recalibration model by interval
    call GATK.ApplyBQSR {
      input:
        input_bam = AddRg.bam_wRG,
        input_bam_index = AddRg.bam_wRG_idx,
        output_bam_basename = sample_name + ".aligned.duplicates_marked.recalibrated",
        recalibration_report = GatherBqsrReports.output_bqsr_report,
        sequence_group_interval = subgroup,
        ref_dict = reference_dict,
        ref_fasta = reference_fasta,
        ref_fasta_index = reference_fasta_index,
        docker_image = gatk_docker,
        gatk_path = gatk_path,
        disk_size = ceil(2.5 * size(AddRg.bam_wRG, "GB")) + 10,
        preemptible_tries = 5
    }
  } 

  # Merge the recalibrated BAM files resulting from by-interval recalibration
  call GATK.GatherBamFiles {
    input:
      input_bams = ApplyBQSR.recalibrated_bam,
      output_bam_basename = sample_name,
      docker_image = gatk_docker,
      gatk_path = gatk_path,
      disk_size = ceil(3 * size(AddRg.bam_wRG, "GB")) + 10,
      preemptible_tries = 5,
      compression_level = 2
  }

  if ( cram_out ) {
    call Bam2Cram {
      input:
        bam = GatherBamFiles.output_bam,
        bai = GatherBamFiles.output_bam_index,
        ref_fa = reference_fasta,
        ref_fai = reference_fasta_index,
        docker = gatk_docker
    }
  }

  # Outputs that will be retained when execution is complete  
  output {
    File duplication_metrics = MarkDuplicates.duplicate_metrics
    File bqsr_report = GatherBqsrReports.output_bqsr_report
    File cram_or_bam = select_first([Bam2Cram.cram, GatherBamFiles.output_bam])
    File cram_or_bam_index = select_first([Bam2Cram.crai, GatherBamFiles.output_bam_index])
  } 
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

# Evenly divide paired fastqs to parallelize alignment
task SplitFastq {
  input {
    File fastq
    Int n_shards
    String docker
  }

  Int disk_gb = ceil(size(fastq, "GB") * 3) + 10
  String out_prefix = basename(fastq, ".fq.gz")

  command <<<
    set -eu -o pipefail

    # Build fastqsplitter command
    cmd="fastqsplitter -i ~{fastq}"
    for i in $( seq 1 ~{n_shards} ); do
      cmd="$cmd -o ~{out_prefix}.$i.fq.gz"
    done

    # Split fastqs
    echo -e "Now splitting fastqs with the following command:\n$cmd\n"
    eval $cmd
  >>>

  output {
    Array[File] fq_shards = glob("~{out_prefix}.*.fq.gz")
  }

  runtime {
    docker: docker
    memory: "3.5GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}

# Align single or paired fastqs with BWA MEM
task Bwa {
  input {
    Array[File] input_fastqs
    String output_basename
    File ref_fasta
    File ref_fasta_index
    File ref_dict
    File? ref_alt
    File ref_amb
    File ref_ann
    File ref_bwt
    File ref_pac
    File ref_sa

    Float mem_gb = 32
    String num_cpu = 8

    Int preemptible_tries = 3
    Int? disk_size

    String docker
  }

  Int default_disk_size = ceil(3 * size(input_fastqs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Align fastq(s) with bwa
    /usr/gitc/bwa mem -v 3 -t ~{num_cpu} -Y \
      ~{ref_fasta} \
      ~{sep=" " input_fastqs} \
    | samtools view -b --output-fmt-option level=2 - \
    > ~{output_basename}.bam
  >>>
  
  output {
    File output_bam = "~{output_basename}.bam"
  }

  runtime {
    preemptible: preemptible_tries
    docker: docker
    memory: "~{mem_gb} GiB"
    cpu: num_cpu
    disks: "local-disk " + select_first([disk_size, default_disk_size]) + " HDD"
  }
}

# Add a read group tag to an aligned BAM
task AddRg {
  input {
    File bam
    File bam_index
    String sample_name
    String gatk_path
    String docker
  }

  String out_filename = basename(bam, ".bam") + ".wRG.bam"
  Int disk_gb = ceil(3 * size(bam, "GB")) + 10

  command {
    set -eu -o pipefail

    ~{gatk_path} \
      AddOrReplaceReadGroups \
      I=~{bam} \
      O=~{out_filename} \
      RGID=1 \
      RGLB=lib1 \
      RGPL=ILLUMINA \
      RGPU=unit1 \
      RGSM=~{sample_name}

    samtools index ~{out_filename}
  }

  output {
    File bam_wRG = "~{out_filename}"
    File bam_wRG_idx = "~{out_filename}.bai"
  }

  runtime {
    preemptible: 5
    docker: docker
    memory: "3.5 GiB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
  }
}

# Sort BAM file by coordinate order and fix tag values for NM and UQ
# Note: taken from GATK-4 WDL but modified to use samtools sort
task SortAndFixTags {
  input {
    File input_bam
    String output_bam_basename
    File ref_dict
    File ref_fasta
    File ref_fasta_index
  
    Int compression_level
    Int preemptible_tries
    Int disk_size
    Int sort_mem_gb = 4
    Int sort_cpu = 4
    Int fix_tags_mem_gb = 2

    String docker_image
    String gatk_path
  }
  Int sort_extra_threads = sort_cpu - 1
  Int total_mem_gb = sort_mem_gb + fix_tags_mem_gb + 2

  command {
    set -o pipefail

    samtools sort \
      -l 1 -O bam -@ ~{sort_extra_threads} \
      ~{input_bam} \
    | ~{gatk_path} \
      --java-options "-Dsamjdk.compression_level=~{compression_level} -Xms~{fix_tags_mem_gb}G" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT ~{output_bam_basename}.bam \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE ~{ref_fasta}
  }
  runtime {
    preemptible: preemptible_tries
    docker: docker_image
    memory: "~{total_mem_gb} GiB"
    cpu: sort_cpu
    disks: "local-disk " + disk_size + " HDD"
  }
  output {
    File output_bam = "~{output_bam_basename}.bam"
    File output_bam_index = "~{output_bam_basename}.bai"
    File output_bam_md5 = "~{output_bam_basename}.bam.md5"
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
      --threads 2 \
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
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 5
  }
}
