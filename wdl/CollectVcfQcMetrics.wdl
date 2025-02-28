# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collect quality metrics for assessing a GATK joint-genotyped VCF (SNVs/indels or SVs)


version 1.0


import "Utilities.wdl" as Utils
import "QcTasks.wdl" as QcTasks


workflow CollectVcfQcMetrics {
  input {
    File vcf
    File vcf_idx

    File? trios_ped_file                      # .ped file of known families for Mendelian transmission analyses
    File? sample_priority_list                # Rank-ordered list of samples to retain for sample-level analyses
    Int n_for_sample_level_analyses = 1000    # Number of samples to use for sample-level summary analyses

    File scatter_intervals_list               # GATK-style intervals file for scattering over vcf (any tabix-compliant interval definitions should work)

    String output_prefix

    String bcftools_docker
    String g2c_analysis_docker
    String linux_docker
  }

  # Get list of samples present in input VCFs
  call Utils.GetSamplesFromVcfHeader as GetSamplesInVcf {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      bcftools_docker = bcftools_docker
  }

  # Clean input .ped file to retain only complete trios
  # TODO: implement this

  # Define list of samples to use for sample-specific analyses
  call ChooseTargetSamples {
    input:
      all_samples_list = GetSamplesInVcf.sample_list,
      sample_priority_list = sample_priority_list,
      n_samples = n_for_sample_level_analyses,
      docker = bcftools_docker
  }

  # Clean scatter intervals
  call Utils.ParseIntervals {
    input:
      intervals_list = scatter_intervals_list,
      docker = linux_docker
  }

  # Parallelize over scatter intervals
  scatter ( interval_info in ParseIntervals.interval_info ) {

    String interval_coords = interval_info.right
    String shard_prefix = basename(vcf, ".vcf.gz") + "." + interval_info.left

    # Extract region, update INFO, and create VCF subsets as needed
    call SliceAndCleanVcf {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        interval = interval_coords,
        target_samples = ChooseTargetSamples.target_samples,
        out_prefix = shard_prefix,
        docker = bcftools_docker
    }

    # Compute site-level metrics
    call QcTasks.CollectSiteMetrics {
      input:
        vcf = SliceAndCleanVcf.sites_vcf,
        vcf_idx = SliceAndCleanVcf.sites_vcf_idx,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Compute site-level benchmarking metrics
    # TODO: implement this

    # Compute sample-level metrics
    # TODO: implement this

    # Compute trio metrics
    # TODO: implement this

    # Compute sample-level benchmarking metrics
  }

  output {}
}


# Define list of "target" samples to use for sample-based analyses
task ChooseTargetSamples {
  input {
    File all_samples_list
    String out_prefix
    File? sample_priority_list
    Int n_samples
    String docker
  }

  String outfile = out_prefix + ".target.samples.list"

  command <<<
    set -eu -o pipefail

    if [ ~{defined(sample_priority_list)} == "true" ]; then
      fgrep \
        -xf ~{all_samples_list} \
        ~{select_first([sample_priority_list])} \
      | head -n ~{n_samples} \
      > "~{outfile}" || true
    else
      shuf \
        --random-source=<( yes 2025 ) \
        ~{all_samples_list} \
      | head -n ~{n_samples} \
      > "~{outfile}" || true
    fi
  >>>

  output {
    File target_samples = outfile
  }

  runtime {
    docker: docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk 25 HDD"
    preemptible: 3
    max_retries: 1
  }
}


# Slice a VCF to a prespecified interval and update all INFO required for QC metric collection
task SliceAndCleanVcf {
  input {
    File vcf
    File vcf_idx
    String interval
    File target_samples
    File? trio_samples
    File? benchmarking_samples
    
    String out_prefix
    
    Int? disk_gb
    Float mem_gb = 4
    Int n_cpu = 2
    
    String docker
  }

  Array[File] all_sample_lists = select_all([trio_samples, target_samples, benchmarking_samples])

  String outfile = out_prefix + ".vcf.gz"
  String sites_outfile = out_prefix + ".sites.vcf.gz"

  Int default_disk_gb = ceil(3 * size(vcf, "GB")) + 10
  Int hdd_gb = select_first([disk_gb, default_disk_gb])

  parameter_meta {
    vcf: {
      localization_optional: true
    }
  }

  command <<<
    set -eu -o pipefail

    # Parse intervals to set bounds on POS (important for avoiding double-counting SVs)
    n_int_parts=$( echo "~{interval}" | sed 's/:\|-/\n/g' | wc -l )
    if [ $n_int_parts -lt 3 ]; then
      min_pos=-1
      max_pos=1000000000
    else
      min_pos=$( echo "~{interval}" | sed 's/:\|-/\t/g' | cut -f2 )
      max_pos=$( echo "~{interval}" | sed 's/:\|-/\t/g' | cut -f3 )
    fi

    # Define superset of samples whose GTs we need for subsequent analysis
    cat ~{sep=" " all_sample_lists} | sort -V | uniq > all.samples.list

    # Symlink vcf_idx to current working dir
    ln -s ~{vcf_idx} .

    # Stream VCF to interval of interest and update INFO
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    bcftools view \
      --regions "~{interval}" \
      ~{vcf} \
    | awk -v min_pos="$min_pos" -v max_pos="$max_pos" \
      '{ if ($1 ~ "^#" || ($2 >= min_pos && $2 <= max_pos)) print }' \
    | bcftools +fill-tags \
    | bcftools view \
      --samples-file all.samples.list \
      --no-update \
    -Oz -o "~{out_prefix}.cleaned_wGTs.vcf.gz"
    tabix -p vcf "~{out_prefix}.cleaned_wGTs.vcf.gz"

    # Make sites VCF for site-level analyses
    bcftools view -G \
      -Oz -o "~{sites_outfile}" \
      "~{out_prefix}.cleaned_wGTs.vcf.gz"
    tabix -p vcf "~{sites_outfile}"

    # Make VCF of only sites & samples appearing in target samples
    # TODO: implement this

    # Make VCF of only sites & samples appearing in complete trios
    # TODO: implement this

    # Make VCF of only sites & samples appearing in benchmarking samples
    # TODO: implement this
  >>>

  output {
    File sites_vcf = sites_outfile
    File sites_vcf_idx = "~{sites_outfile}.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + hdd_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}

