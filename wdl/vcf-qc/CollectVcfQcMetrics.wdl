# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collect quality metrics for assessing a GATK joint-genotyped VCF (SNVs/indels and/or SVs)


version 1.0


import "Utilities.wdl" as Utils
import "QcTasks.wdl" as QcTasks


workflow CollectVcfQcMetrics {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    File? trios_fam_file                       # .fam file of known families for Mendelian transmission analyses
    File? sample_priority_list                 # Rank-ordered list of samples to retain for sample-level analyses
    Int n_for_sample_level_analyses = 1000     # Number of samples to use for sample-level summary analyses

    Boolean shard_vcf = true                   # Should the input VCF be sharded for QC collection?
    File? scatter_intervals_list               # GATK-style intervals file for scattering over vcf 
                                               # (any tabix-compliant interval definitions should work)
    Int n_records_per_shard = 25000            # Number of records per shard. This will only be used as a backup if 
                                               # scatter_intervals_list is not provided and shard_vcf is true

    Float common_af_cutoff = 0.001             # Minimum AF for a variant to be included in common variant subsets

    String output_prefix

    String bcftools_docker
    String g2c_analysis_docker
    String linux_docker
  }

  #####################
  ### SAMPLE MANAGEMENT
  #####################

  # Get list of samples present in input VCFs
  call Utils.GetSamplesFromVcfHeader as GetSamplesInVcf {
    input:
      vcf = vcfs[0],
      vcf_idx = vcf_idxs[0],
      bcftools_docker = bcftools_docker
  }

  # Clean input .fam file to retain only complete trios
  if( defined(trios_fam_file) ) {
    call CleanFam {
      input:
        fam_file = select_first([trios_fam_file]),
        all_samples_list = GetSamplesInVcf.sample_list,
        docker = g2c_analysis_docker
    }
  }
  File all_qc_samples_list = select_first([CleanFam.all_samples_no_probands_list,
                                           GetSamplesInVcf.sample_list])

  # Define list of samples to use for sample-specific analyses
  call ChooseTargetSamples {
    input:
      all_samples_list = all_qc_samples_list,
      sample_priority_list = sample_priority_list,
      n_samples = n_for_sample_level_analyses,
      out_prefix = output_prefix,
      docker = bcftools_docker
  }

  #####################
  ### VCF PREPROCESSING
  #####################

  # Read scatter intervals, if optioned
  if ( defined(scatter_intervals_list) ) {
    call Utils.ParseIntervals {
      input:
        intervals_list = select_first([scatter_intervals_list]),
        docker = linux_docker
    }
  }

  # Preprocess each VCF according to desired scatter behavior
  Array[Pair[File, File]] input_vcf_infos = zip(vcfs, vcf_idxs)
  scatter ( input_vcf_info in input_vcf_infos ) {

    File vcf = input_vcf_info.left
    File vcf_idx = input_vcf_info.right

    if ( shard_vcf ) {

      # Default to scattering over prespecified intervals, as this is most compute efficient
      if ( defined(scatter_intervals_list) ) {        
        scatter ( interval_info in select_first([ParseIntervals.interval_info]) ) {

          String interval_coords = interval_info[1]
          String shard_prefix = basename(vcf, ".vcf.gz") + "." + interval_info[0]

          # Extract desired interval from VCF
          call Utils.StreamSliceVcf as SliceInterval {
            input:
              vcf = vcf,
              vcf_idx = vcf_idx, 
              interval = interval_coords,
              outfile_name = shard_prefix + ".vcf.gz"
          }
        }
      }

      # If scatter_intervals_list is not provided, shard the VCF the old fashioned way
      if ( !defined(scatter_intervals_list) ) {

        # Shard VCF on a VM (costly & slow for large VCFs)
        call Utils.ShardVcf {
          input:
            vcf = vcf,
            vcf_idx = vcf_idx,
            records_per_shard = n_records_per_shard,
            bcftools_docker = bcftools_docker,
            n_preemptible = 1
        }
      }
    }

    # Unify the three approaches for sharding a VCF prior to preprocessing
    Array[File] vcf_shards = select_first([SliceInterval.vcf_slice,
                                           ShardVcf.vcf_shards,
                                           [vcf]])
    Array[File] vcf_shard_idxs = select_first([SliceInterval.vcf_slice_idx,
                                               ShardVcf.vcf_shard_idxs,
                                               [vcf_idx]])


    # Scatter over VCF shards and preprocess each shard
    Array[Pair[File, File]] pp_vcf_infos = zip(vcf_shards, vcf_shard_idxs)
    scatter ( pp_vcf_info in pp_vcf_infos ) {

      File pp_vcf = pp_vcf_info.left
      File pp_vcf_idx = pp_vcf_info.right

      call PreprocessVcf {
        input:
          vcf = vcf,
          vcf_idx = vcf_idx,
          target_samples = ChooseTargetSamples.target_samples,
          out_prefix = basename(vcf, ".vcf.gz"),
          docker = bcftools_docker
      }
    }
  }

  #####################
  ### METRIC COLLECTION
  #####################

  Array[File] site_vcf_shards = flatten(PreprocessVcf.sites_vcf)
  Array[File] site_vcf_shard_idxs = flatten(PreprocessVcf.sites_vcf_idx)
  Array[Pair[File, File]] site_vcf_info = zip(site_vcf_shards, site_vcf_shard_idxs)

  # Collect site-level metrics for all preprocessed shards
  scatter ( shard_info in site_vcf_info ) {

    # Compute site-level metrics
    call QcTasks.CollectSiteMetrics {
      input:
        vcf = shard_info.left,
        vcf_idx = shard_info.right,
        n_samples = ChooseTargetSamples.n_samples_all,
        common_af_cutoff = common_af_cutoff,
        g2c_analysis_docker = g2c_analysis_docker
    }

    # Compute site-level benchmarking metrics
    # TODO: implement this

    # Compute sample-level metrics
    # TODO: implement this

    # Compute trio metrics
    # TODO: implement this

    # Compute sample-level benchmarking metrics
    # TODO: implement this
  }

  ##################
  ### OUTPUT CLEANUP
  ##################

  # Collapse SV sites
  # TODO: implement this

  # Collapse size distributions
  call QcTasks.SumCompressedDistribs as SumSizeDistribs {
    input:
      distrib_tsvs = CollectSiteMetrics.size_distrib,
      out_prefix = output_prefix + ".size_distribution",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collapse freq distributions
  call QcTasks.SumCompressedDistribs as SumAfDistribs {
    input:
      distrib_tsvs = CollectSiteMetrics.af_distrib,
      out_prefix = output_prefix + ".af_distribution",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collapse size X freq distributions
  call QcTasks.SumCompressedDistribs as SumJointDistribs {
    input:
      distrib_tsvs = CollectSiteMetrics.size_vs_af_distrib,
      n_key_columns = 3,
      out_prefix = output_prefix + ".size_vs_af_distribution",
      g2c_analysis_docker = g2c_analysis_docker
  }

  output {
    # For now, just outputting the merged distribution files
    # This will eventually be extended
    File size_distrib = SumSizeDistribs.merged_distrib
    File af_distrib = SumAfDistribs.merged_distrib
    File size_vs_af_distrib = SumJointDistribs.merged_distrib
  }
}


# Clean input .fam file to restrict to complete trios present in input VCF
task CleanFam {
  input {
    File fam_file
    File all_samples_list
    String docker
  }

  String fam_out_fname = basename(fam_file, ".fam") + ".cleaned_trios.fam"
  String proband_out_fname = basename(fam_file, ".fam") + ".proband_ids.list"
  String unrelated_out_fname = basename(all_samples_list) + ".no_probands.list"
  
  command <<<
    set -eu -o pipefail

    # Clean .fam file to duos or trios present in all_samples_list
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/subset_fam.R \
      --in-fam ~{fam_file} \
      --all-samples ~{all_samples_list} \
      --out-fam duos_and_trios.fam

    # Write list of all probands (can exclude these samples for better AF estimates)
    cut -f1 duos_and_trios.fam > "~{proband_out_fname}"

    # Further restrict .fam to complete trios
    awk -v FS="\t" -v OFS="\t" \
      '{ if ($3!=0 && $4!=0) print }' \
      duos_and_trios.fam \
    > "~{fam_out_fname}"

    # Remove probands from all samples
    fgrep \
      -xvf "~{proband_out_fname}" \
      ~{all_samples_list} \
    > "~{unrelated_out_fname}"
  >>>

  output {
    File trios_fam = "~{fam_out_fname}"
    File probands_list = "~{proband_out_fname}"
    File all_samples_no_probands_list = "~{unrelated_out_fname}"
  }

  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk 20 HDD"
    preemptible: 3
    max_retries: 1
  }
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

    touch "~{outfile}"

    # First, add samples from a priority list, if optioned
    if [ ~{defined(sample_priority_list)} == "true" ]; then
      fgrep \
        -xf ~{all_samples_list} \
        ~{select_first([sample_priority_list])} \
      | head -n ~{n_samples} \
      >> "~{outfile}" || true
      fgrep \
        -xvf ~{select_first([sample_priority_list])} \
        ~{all_samples_list} \
      > remainder.samples.list
    else
      cp ~{all_samples_list} remainder.samples.list
    fi

    # Second, supplement from random selection of remaining samples
    subtotal=$( cat ~{outfile} | wc -l )
    n_to_add=$(( ~{n_samples} - $subtotal ))
    if [ $n_to_add -gt 0 ]; then
      shuf \
        --random-source=<( yes 2025 ) \
        remainder.samples.list \
      | head -n $n_to_add \
      >> "~{outfile}" || true
    fi
  >>>

  output {
    File target_samples = outfile
    Int n_samples_all = length(read_lines(all_samples_list))
  }

  runtime {
    docker: docker
    memory: "3.75 GB"
    cpu: 2
    disks: "local-disk 20 HDD"
    preemptible: 3
    max_retries: 1
  }
}


# Update all VCF INFO required for QC metric collection, make all necessary genotype subsets
task PreprocessVcf {
  input {
    File vcf
    File vcf_idx
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

  command <<<
    set -eu -o pipefail

    # Define superset of samples whose GTs we need for subsequent analysis
    cat ~{sep=" " all_sample_lists} | sort -V | uniq > all.samples.list

    # Preprocess VCF
    bcftools +fill-tags ~{vcf} -- -t AN,AC,AF,AC_Hemi,AC_Het,AC_Hom,HWE \
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

