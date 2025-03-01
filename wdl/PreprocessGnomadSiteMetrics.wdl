# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collect variant summary metrics for gnomAD to match G2C QC
# Designed to run on a single chromosome (based on how gnomAD data are stored online)


version 1.0


import "Utilities.wdl" as Utils
import "QcTasks.wdl" as QcTasks


workflow PreprocessGnomadSiteMetrics {
  input {
    File snv_vcf
    File snv_vcf_idx
    File snv_scatter_intervals
    Int snv_n_samples

    File sv_vcf
    File sv_vcf_idx
    File sv_scatter_intervals
    Int sv_n_samples

    String output_prefix

    String linux_docker
    String bcftools_docker
    String g2c_analysis_docker
  }

  # Clean SNV scatter intervals
  call Utils.ParseIntervals as MakeSnvIntervals {
    input:
      intervals_list = snv_scatter_intervals,
      docker = linux_docker
  }

  # Scatter SNV VCF over intervals
  scatter ( snv_interval_info in MakeSnvIntervals.interval_info ) {

    String snv_interval_coords = snv_interval_info.right
    String snv_shard_prefix = output_prefix + ".snv." + snv_interval_info.left

    # Slice gnomAD SNV VCF to interval and drop unnecessary annotations
    call SliceVcf as SliceSnvShard {
      input:
        vcf = snv_vcf,
        vcf_idx = snv_vcf_idx,
        interval = snv_interval_coords,
        out_prefix = snv_shard_prefix,
        bcftools_docker = bcftools_docker
    }

    # Collect site metrics for SNVs
    call QcTasks.CollectSiteMetrics as CollectSnvMetrics {
      input:
        vcf = SliceSnvShard.cleaned_vcf,
        vcf_idx = SliceSnvShard.cleaned_vcf_idx,
        n_samples = snv_n_samples,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # Concatenate sites BEDs for SNVs
  call Utils.ConcatTextFiles as ConcatSnvs {
    input:
      shards = select_all(CollectSnvMetrics.snv_sites),
      concat_command = "zcat",
      sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
      compression_command = "bgzip -c",
      input_has_header = true,
      output_filename = output_prefix + ".snv.sites.bed.gz",
      docker = bcftools_docker
  }

  # Concatenate sites BEDs for indels
  call Utils.ConcatTextFiles as ConcatIndels {
    input:
      shards = select_all(CollectSnvMetrics.indel_sites),
      concat_command = "zcat",
      sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
      compression_command = "bgzip -c",
      input_has_header = true,
      output_filename = output_prefix + ".indel.sites.bed.gz",
      docker = bcftools_docker
  }

  # Clean SV scatter intervals
  call Utils.ParseIntervals as MakeSvIntervals {
    input:
      intervals_list = sv_scatter_intervals,
      docker = linux_docker
  }

  # Scatter SV VCF over intervals
  scatter ( sv_interval_info in MakeSvIntervals.interval_info ) {

    String sv_interval_coords = sv_interval_info.right
    String sv_shard_prefix = output_prefix + ".snv." + sv_interval_info.left

    # Slice gnomAD SV VCF to interval and drop unnecessary annotations
    call SliceVcf as SliceSvShard {
      input:
        vcf = sv_vcf,
        vcf_idx = sv_vcf_idx,
        interval = sv_interval_coords,
        out_prefix = sv_shard_prefix,
        bcftools_docker = bcftools_docker
    }

    # Collect site metrics for SVs
    call QcTasks.CollectSiteMetrics as CollectSvMetrics {
      input:
        vcf = SliceSvShard.cleaned_vcf,
        vcf_idx = SliceSvShard.cleaned_vcf_idx,
        n_samples = snv_n_samples,
        g2c_analysis_docker = g2c_analysis_docker
    }
  }

  # Concatenate sites BEDs for SVs
  call Utils.ConcatTextFiles as ConcatSvs {
    input:
      shards = select_all(CollectSnvMetrics.sv_sites),
      concat_command = "zcat",
      sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
      compression_command = "bgzip -c",
      input_has_header = true,
      output_filename = output_prefix + ".sv.sites.bed.gz",
      docker = bcftools_docker
  }

  # Collapse size distributions
  call QcTasks.SumCompressedDistribs as SumSizeDistribs {
    input:
      distrib_tsvs = flatten([CollectSnvMetrics.size_distrib, 
                              CollectSvMetrics.size_distrib]),
      out_prefix = output_prefix + ".size_distribution",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collapse AF distributions
  call QcTasks.SumCompressedDistribs as SumAfDistribs {
    input:
      distrib_tsvs = flatten([CollectSnvMetrics.af_distrib, 
                              CollectSvMetrics.af_distrib]),
      out_prefix = output_prefix + ".af_distribution",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collapse 2D size vs. AF distributions
  call QcTasks.SumCompressedDistribs as SumJointDistribs {
    input:
      distrib_tsvs = flatten([CollectSnvMetrics.size_vs_af_distrib, 
                              CollectSvMetrics.size_vs_af_distrib]),
      n_key_columns = 3,
      out_prefix = output_prefix + ".size_vs_af_distribution",
      g2c_analysis_docker = g2c_analysis_docker
  }

  output {
    File snv_sites_bed = ConcatSnvs.merged_file
    File indel_sites_bed = ConcatIndels.merged_file
    File sv_sites_bed = ConcatSvs.merged_file
    File size_distrib = SumSizeDistribs.merged_distrib
    File af_distrib = SumAfDistribs.merged_distrib
    File size_vs_af_distrib = SumJointDistribs.merged_distrib
  }
}


# Slice & clean a gnomAD VCF to minimal info for site metric collection
task SliceVcf {
  input {
    File vcf
    File vcf_idx
    String interval

    String out_prefix

    Int disk_gb = 25
    String bcftools_docker
  }

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

    # Symlink vcf_idx to current working dir
    ln -s ~{vcf_idx} .

    # Stream VCF to interval of interest and remove all INFO except for 
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    bcftools view \
      --regions "~{interval}" \
      ~{vcf} \
    | awk -v min_pos="$min_pos" -v max_pos="$max_pos" \
      '{ if ($1 ~ "^#" || ($2 >= min_pos && $2 <= max_pos)) print }' \
    | bcftools annotate \
      -x "^INFO/END,INFO/SVTYPE,INFO/SVLEN,INFO/AC,INFO/AF,INFO/HWE,INFO/ExcHet" \
      -Oz -o "~{out_prefix}.cleaned.vcf.gz"
    tabix -p vcf "~{out_prefix}.cleaned.vcf.gz"

  >>>

  output {
    File cleaned_vcf = "~{out_prefix}.cleaned.vcf.gz"
    File cleaned_vcf_idx = "~{out_prefix}.cleaned.vcf.gz.tbi"
  }

  runtime {
    docker: bcftools_docker
    memory: "2 GB"
    cpu: 1
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}

