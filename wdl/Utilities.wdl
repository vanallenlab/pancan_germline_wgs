# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Common WDL tasks shared across workflows


version 1.0


# Compute a BEDGraph of feature density from one or more BED files
task CalcBedDensity {
  input {
    Array[File] beds
    String bed_concat_cmd = "zcat"

    File genome_file # BEDTools-style genome file

    String output_prefix
    
    String bedtools_docker
  }

  Int disk_gb = ceil(2 * size(beds, "GB")) + 10
  File outfile = output_prefix + ".density.bed.gz"

  command <<<
    set -eu -o pipefail

    cat ~{write_lines(beds)} \
    | ~{bed_concat_cmd} - \
    | awk -v OFS="\t" '{ print $1, $2, $3 }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools genomecov -bg -g ~{genome_file} -i - \
    | bgzip -c \
    > "~{outfile}"
  >>>

  output {
    File density_bed = "~{outfile}"
  }

  runtime {
    docker: bedtools_docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 1
    maxRetries: 1
  }
}


task ConcatVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs
    String out_prefix

    String bcftools_concat_options = ""

    Float mem_gb = 3.5
    Int cpu_cores = 2
    Int? disk_gb

    String bcftools_docker
  }

  String out_filename = out_prefix + ".vcf.gz"

  Int default_disk_gb = ceil(2.5 * size(vcfs, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools concat \
      ~{bcftools_concat_options} \
      --file-list ~{write_lines(vcfs)} \
      -O z \
      -o ~{out_filename} \
      --threads ~{cpu_cores}

    tabix -p vcf -f ~{out_filename}
  >>>

  output {
    File merged_vcf = "~{out_filename}"
    File merged_vcf_idx = "~{out_filename}.tbi"
  }

  runtime {
    docker: bcftools_docker
    memory: mem_gb + " GB"
    cpu: cpu_cores
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    preemptible: 3
  }
}


task CopyGcpObjects {
  input {
    Array[File] files_to_copy
    String destination
    String gsutil_cp_options = ""
  }

  Int disk_gb = ceil(1.2 * size(files_to_copy, "GB")) + 10

  command <<<
    set -eu -o pipefail

    gsutil -m cp ~{gsutil_cp_options} ~{sep=" " files_to_copy} "~{destination}"
  >>>

  output {}

  runtime {
    docker: "google/cloud-sdk"
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task CountRecordsInVcf {
  input {
    File vcf
    File? vcf_idx
    String bcftools_docker
  }

  Int disk_gb = ceil(1.2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    if [ ~{defined(vcf_idx)} == "false" ]; then
      tabix -p vcf -f ~{vcf}
    fi

    bcftools query -f '%CHROM\n' ~{vcf} | wc -l > record_count.txt
  >>>

  output {
    Int n_records = read_int("record_count.txt")
  }

  runtime {
    docker: bcftools_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


# Generic task to efficiently download a file from an FTP server
task FtpDownload {
  input {
    String ftp_url
    String output_name

    String lftp_docker
    
    Int max_download_tries = 3
    Int disk_gb = 150
    Int n_cpu = 4
    Float mem_gb = 8
  }

  command {
    set -euo pipefail

    echo "Downloading ${ftp_url} to ${output_name}"
    
    # Retry loop with up to 5 attempts
    count=0
    success=0

    while [ $count -lt ~{max_download_tries} ]; do
      echo "Attempt $((count+1))..."
      lftp -c "set net:max-retries 2; set net:timeout 30; pget -n 8 -c ${ftp_url} -o ${output_name}" && success=1 && break
      count=$((count+1))
      echo "Attempt $count failed. Retrying..."
      sleep 30
    done

    if [ $success -ne 1 ]; then
      echo "Failed to download after ~{max_download_tries} attempts."
      exit 1
    fi

    echo "Download completed successfully."
  }

  output {
    File downloaded_file = output_name
  }

  runtime {
    docker: lftp_docker
    disks: "local-disk " + disk_gb + " HDD"
    cpu: n_cpu
    memory: mem_gb + " GB"
    preemptible: 0
  }
}


task GetContigsFromFai {
  input {
    File ref_fai
    String docker
  }

  command <<<
    set -eu -o pipefail

    cut -f1 ~{ref_fai} > contigs.list
  >>>

  output {
    Array[String] contigs = read_lines("contigs.list")
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 10 HDD"
    preemptible: 3
  }
}


task GetContigsFromVcfHeader {
  input {
    File vcf
    File vcf_idx
    String docker
  }

  Int disk_gb = ceil(1.2 * size(vcf, "GB")) + 10

  parameter_meta {
    vcf: {
      localization_optional: true
    }
  }

  command <<<
    set -eu -o pipefail

    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    ln -s ~{vcf_idx} .

    tabix -H ~{vcf} \
    | fgrep "##contig" \
    | sed 's/ID=/\t/g' \
    | cut -f2 \
    | cut -f1 -d, \
    | sort -V \
    > contigs.list
  >>>

  output {
    Array[String] contigs = read_lines("contigs.list")
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task IntersectTextFiles {
  input {
    Array[File] files
    String outfile = "intersection.txt"
    String docker
  }

  Int disk_gb = 2 * ceil(size(files, "GB")) + 10

  command <<<
    set -eu -o pipefail

    cat ~{files[0]} > "~{outfile}"
    while read fid; do
      fgrep -xf $fid "~{outfile}" > "~{outfile}2"
      mv "~{outfile}2" "~{outfile}"
    done < ~{write_lines(files)}
  >>>

  output {
    File intersection_file = "~{outfile}"
  }

  runtime {
    docker: docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task MakeTabixIndex {
  input {
    File input_file
    String file_type = "vcf"
    String docker
    Int n_cpu = 2
    Float mem_gb = 3.5
  }

  String outfile = basename(input_file, "gz") + "gz.tbi"
  Int disk_gb = ceil(1.25 * size(input_file, "GB")) + 20

  command <<<
    set -eu -o pipefail

    # Due to weirdness with call cacheing and where tabix writes 
    # indexes by default, it seems safest to relocate the input
    # file to the pwd before indexing it
    
    mv ~{input_file} ./
    tabix -p ~{file_type} -f ~{basename(input_file)}
  >>>

  output {
    File tbi = "~{outfile}"
  }

  runtime {
    docker: docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


task ShardVcf {
  input {
    File vcf
    File vcf_idx
    Int records_per_shard
    String bcftools_docker
    Int? disk_gb
    Int n_preemptible = 3
  }

  String out_prefix = basename(vcf, ".vcf.gz") + ".sharded"
  Int use_disk_gb = select_first([disk_gb, ceil(5 * size(vcf, "GB")) + 20])

  command <<<
    set -eu -o pipefail

    # Make an empty shard in case the input VCF is totally empty
    bcftools view -h ~{vcf} | bgzip -c > "~{out_prefix}.0.vcf.gz"

    bcftools +scatter \
      -O z3 -o . -p "~{out_prefix}". \
      -n ~{records_per_shard} \
      ~{vcf}

    # Print all VCFs to stdout for logging purposes
    find ./ -name "*.vcf.gz"

    # Index all shards
    find ./ -name "~{out_prefix}.*.vcf.gz" \
    | xargs -I {} tabix -p vcf -f {}
  >>>

  output {
    Array[File] vcf_shards = glob("~{out_prefix}.*.vcf.gz")
    Array[File] vcf_shard_idxs = glob("~{out_prefix}.*.vcf.gz.tbi")
  }

  runtime {
    cpu: 2
    memory: "3.75 GiB"
    disks: "local-disk " + use_disk_gb + " HDD"
    bootDiskSizeGb: 10
    docker: bcftools_docker
    preemptible: n_preemptible
    maxRetries: 1
  }
}


# Atomizes a GATK-style interval list into single-interval files
task SplitIntervalList {
  input {
    File interval_list
    String linux_docker = "marketplace.gcr.io/google/ubuntu1804"
    Int n_preemptible = 3
  }

  command <<<
    set -eu -o pipefail

    mkdir scatterDir

    fgrep "@" ~{interval_list} > header.txt || true

    split -l 1 -a 6 -d <( fgrep -v "@" ~{interval_list} ) scatterDir/shard

    while read shard; do
      i=$( basename $shard | sed 's/^shard//g' | awk '{ printf "%06d\n", $1+1 }' )
      cat header.txt $shard > scatterDir/$i-scattered.interval_list
      rm $shard
    done < <( find scatterDir/ -name "shard*" | sort -n )
  >>>

  output {
    Array[File] output_intervals = glob("scatterDir/*-scattered.interval_list")
  }

  runtime {
    cpu: 1
    memory: "1.75 GiB"
    disks: "local-disk 25 HDD"
    docker: linux_docker
    preemptible: n_preemptible
    maxRetries: 1
  }
}


# Checks a VCF header for fields indicating that it might contain mCNVs
# Functionally equivalent to McnvHeaderCheck, but uses htslib remote streaming
# Unable to disable localization_optional in parameter_meta so we must 
# have two versions of roughly the same function
task StreamedMcnvHeaderCheck {
  input {
    File vcf
    File vcf_idx
    String docker
  }

  parameter_meta {
    vcf: {
      localization_optional: true
    }
  }

  command <<<
    set -eu -o pipefail

    ln -s ~{vcf_idx} .
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    bcftools view --header-only ~{vcf} > header.vcf

    # Check for header rows potentially indicative of MCNVs
    touch mcnv.header.vcf
    fgrep "ALT=<ID=CNV" header.vcf > mcnv.header.vcf || true
    fgrep "ALT=<ID=MCNV" header.vcf >> mcnv.header.vcf || true
    fgrep "FILTER=<ID=MULTIALLELIC" header.vcf >> mcnv.header.vcf || true
    if [ $( cat mcnv.header.vcf | wc -l ) -gt 0 ]; then
      echo "true" > has_mcnvs.txt
    else
      echo "false" > has_mcnvs.txt
    fi
  >>>

  output {
    Boolean has_mcnvs = read_boolean("has_mcnvs.txt")
  }

  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk 15 HDD"
    preemptible: 3
    maxRetries: 2
  }
}


task StreamSamplesFromVcfHeader {
  input {
    File vcf
    File vcf_idx
    String bcftools_docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + ".samples.list"
  Int disk_gb = ceil(1.2 * size(vcf, "GB")) + 10

  parameter_meta {
    vcf: {
      localization_optional: true
    }
  }

  command <<<
    set -eu -o pipefail

    ln -s ~{vcf_idx} .
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    bcftools query -l ~{vcf} > ~{out_filename}
  >>>

  output {
    String sample_list = out_filename
    Int n_samples = length(read_lines(out_filename))
  }

  runtime {
    docker: bcftools_docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}


task SumSvCountsPerSample {
  input {
    Array[File] count_tsvs   # Expects svtk count-svtypes output format
    String output_prefix

    String docker
    Float mem_gb = 3.75
    Int n_cpu = 2
  }

  Int disk_gb = ceil(2 * size(count_tsvs, "GB")) + 10
  String outfile = output_prefix + ".counts.tsv"

  command <<<
    set -euo pipefail

    /opt/pancan_germline_wgs/scripts/gatksv_helpers/sum_svcounts.py \
      --outfile "~{outfile}" \
      ~{sep=" " count_tsvs}
  >>>

  output {
    File summed_tsv = "~{outfile}"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
  }
}
