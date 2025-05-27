# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform two-sided site-level benchmarking of a "source" callset versus a reference "target" callset
# Note: this is a subroutine called within BenchmarkSites.wdl, and is designed to be
# invoked by this parent workflow

# Expects that both callsets have already been processed by CollectVcfQcMetrics.wdl


version 1.0


import "QcTasks.wdl" as QcTasks
import "Utilities.wdl" as Utils


workflow BenchmarkSitesSingle {
  input {
    File eval_interval_bed
    File genome_file

    File? source_snv_bed
    File? source_snv_bed_idx
    File? target_snv_bed
    File? target_snv_bed_idx

    File? source_indel_bed
    File? source_indel_bed_idx
    File? target_indel_bed
    File? target_indel_bed_idx

    File? source_sv_bed
    File? source_sv_bed_idx
    File? target_sv_bed
    File? target_sv_bed_idx

    String eval_prefix
    String source_prefix
    String target_prefix

    Float common_af_cutoff

    Int total_shards

    String bcftools_docker
    String g2c_analysis_docker
  }

  String ppv_prefix = eval_prefix + "." + source_prefix + "_vs_" + target_prefix
  String sens_prefix = eval_prefix + "." + target_prefix + "_vs_" + source_prefix

  Boolean do_snv = if defined(source_snv_bed) && defined(source_snv_bed_idx)
                   && defined(target_snv_bed) && defined(target_snv_bed_idx)
                   then true else false

  Boolean do_indel = if defined(source_indel_bed) && defined(source_indel_bed_idx)
                     && defined(target_indel_bed) && defined(target_indel_bed_idx)
                     then true else false
  Boolean do_sv = if defined(source_sv_bed) && defined(source_sv_bed_idx)
                  && defined(target_sv_bed) && defined(target_sv_bed_idx) 
                  then true else false
  Int n_var_types = (if do_snv then 1 else 0) + 
                    (if do_indel then 1 else 0) + 
                    (if do_sv then 1 else 0)

  # Shard evaluation intervals, balancing on genomic bp
  Int target_shards = floor(total_shards / (2 * n_var_types))
  call QcTasks.ShardIntervals {
    input:
      intervals_bed = eval_interval_bed,
      n_shards = target_shards,
      prefix = eval_prefix,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Process each eval interval in parallel
  scatter ( eval_shard in ShardIntervals.interval_shards ) {

    String shard_prefix = basename(eval_shard, ".bed.gz")

    # Compare SNVs, if provided
    if ( do_snv ) {
      call PrepSites as PrepSourceSnvs {
        input:
          beds = select_all([source_snv_bed]),
          bed_idxs = select_all([source_snv_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 0,
          max_size = 0,
          prefix = source_prefix + ".snv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call PrepSites as PrepTargetSnvs {
        input:
          beds = select_all([target_snv_bed]),
          bed_idxs = select_all([target_snv_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 0,
          max_size = 0,
          prefix = target_prefix + ".snv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareSnvPpv {
        input:
          query_bed = PrepSourceSnvs.query_bed,
          ref_bed = PrepTargetSnvs.ref_bed,
          genome_file = genome_file,
          mode = "exact",
          common_af = common_af_cutoff,
          prefix = ppv_prefix + ".snv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareSnvSens {
        input:
          query_bed = PrepTargetSnvs.query_bed,
          ref_bed = PrepSourceSnvs.ref_bed,
          genome_file = genome_file,
          mode = "exact",
          common_af = common_af_cutoff,
          degenerate = true,
          prefix = sens_prefix + ".snv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }

    # Compare indels, if provided
    if ( do_indel ) {
      call PrepSites as PrepSourceIndels {
        input:
          beds = select_all([source_indel_bed, source_sv_bed]),
          bed_idxs = select_all([source_indel_bed_idx, source_sv_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 1,
          max_size = 49,
          prefix = source_prefix + ".indel." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call PrepSites as PrepTargetIndels {
        input:
          beds = select_all([target_indel_bed, target_sv_bed]),
          bed_idxs = select_all([target_indel_bed_idx, target_sv_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 1,
          max_size = 49,
          prefix = target_prefix + ".indel." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareIndelPpv {
        input:
          query_bed = PrepSourceIndels.query_bed,
          ref_bed = PrepTargetIndels.ref_bed,
          genome_file = genome_file,
          mode = "both",
          common_af = common_af_cutoff,
          prefix = ppv_prefix + ".indel." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareIndelSens {
        input:
          query_bed = PrepTargetIndels.query_bed,
          ref_bed = PrepSourceIndels.ref_bed,
          genome_file = genome_file,
          mode = "both",
          common_af = common_af_cutoff,
          degenerate = true,
          prefix = sens_prefix + ".indel." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }

    # Compare SVs, if provided
    if ( do_sv ) {
      call PrepSites as PrepSourceSvs {
        input:
          beds = select_all([source_sv_bed, source_indel_bed]),
          bed_idxs = select_all([source_sv_bed_idx, source_indel_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 50,
          max_size = 1000000000,
          prefix = source_prefix + ".sv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call PrepSites as PrepTargetSvs {
        input:
          beds = select_all([target_sv_bed, target_indel_bed]),
          bed_idxs = select_all([target_sv_bed_idx, target_indel_bed_idx]),
          eval_interval_bed = eval_shard,
          min_size = 50,
          max_size = 1000000000,
          prefix = target_prefix + ".sv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareSvPpv {
        input:
          query_bed = PrepSourceSvs.query_bed,
          ref_bed = PrepTargetSvs.ref_bed,
          genome_file = genome_file,
          mode = "both",
          overlap_pad = 2,
          common_af = common_af_cutoff,
          prefix = ppv_prefix + ".sv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }

      call CompareSites as CompareSvSens {
        input:
          query_bed = PrepTargetSvs.query_bed,
          ref_bed = PrepSourceSvs.ref_bed,
          genome_file = genome_file,
          mode = "both",
          overlap_pad = 2,
          common_af = common_af_cutoff,
          degenerate = true,
          prefix = sens_prefix + ".sv." + shard_prefix,
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Collapse SNV common site comparisons
  Array[File] common_snv_ppv_beds = select_all(CompareSnvPpv.common_sites_bed)
  if ( length(common_snv_ppv_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseCommonSnvPpvSites {
      input:
        shards = common_snv_ppv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = ppv_prefix + ".common_sites.snvs.bed.gz",
        docker = bcftools_docker
      }
  }
  Array[File] common_snv_sens_beds = select_all(CompareSnvSens.common_sites_bed)
  if ( length(common_snv_sens_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseCommonSnvSensSites {
      input:
        shards = common_snv_sens_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = sens_prefix + ".common_sites.snvs.bed.gz",
        docker = bcftools_docker
      }
  }

  # Collapse indel common site comparisons
  Array[File] common_indel_ppv_beds = select_all(CompareIndelPpv.common_sites_bed)
  if ( length(common_indel_ppv_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseCommonIndelPpvSites {
      input:
        shards = common_indel_ppv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = ppv_prefix + ".common_sites.indels.bed.gz",
        docker = bcftools_docker
      }
  }
  Array[File] common_indel_sens_beds = select_all(CompareIndelSens.common_sites_bed)
  if ( length(common_indel_sens_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseCommonIndelSensSites {
      input:
        shards = common_indel_sens_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = sens_prefix + ".common_sites.indels.bed.gz",
        docker = bcftools_docker
      }
  }

  # Collapse SV common site comparisons
  Array[File] common_sv_ppv_beds = select_all(CompareSvPpv.common_sites_bed)
  if ( length(common_sv_ppv_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseCommonSvPpvSites {
      input:
        shards = common_sv_ppv_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = ppv_prefix + ".common_sites.svs.bed.gz",
        docker = bcftools_docker
      }
  }
  Array[File] common_sv_sens_beds = select_all(CompareSvSens.common_sites_bed)
  if ( length(common_sv_sens_beds) > 0 ) {
    call Utils.ConcatTextFiles as CollapseCommonSvSensSites {
      input:
        shards = common_sv_sens_beds,
        concat_command = "zcat",
        sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
        compression_command = "bgzip -c",
        input_has_header = true,
        output_filename = sens_prefix + ".common_sites.svs.bed.gz",
        docker = bcftools_docker
      }
  }

  # Sum compressed concordance results by size
  call QcTasks.SumCompressedDistribs as SumPpvBySize {
    input:
      distrib_tsvs = flatten([select_all(CompareSnvPpv.val_by_size_tsv),
                              select_all(CompareIndelPpv.val_by_size_tsv),
                              select_all(CompareSvPpv.val_by_size_tsv)]),
      n_key_columns = 3,
      out_prefix = ppv_prefix + ".concordance_by_size",
      g2c_analysis_docker = g2c_analysis_docker
  }
  call QcTasks.SumCompressedDistribs as SumSensBySize {
    input:
      distrib_tsvs = flatten([select_all(CompareSnvSens.val_by_size_tsv),
                              select_all(CompareIndelSens.val_by_size_tsv),
                              select_all(CompareSvSens.val_by_size_tsv)]),
      n_key_columns = 3,
      out_prefix = sens_prefix + ".concordance_by_size",
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Sum compressed concordance results by frequency
  call QcTasks.SumCompressedDistribs as SumPpvByFreq {
    input:
      distrib_tsvs = flatten([select_all(CompareSnvPpv.val_by_freq_tsv),
                              select_all(CompareIndelPpv.val_by_freq_tsv),
                              select_all(CompareSvPpv.val_by_freq_tsv)]),
      n_key_columns = 3,
      out_prefix = ppv_prefix + ".concordance_by_af",
      g2c_analysis_docker = g2c_analysis_docker
  }
  call QcTasks.SumCompressedDistribs as SumSensByFreq {
    input:
      distrib_tsvs = flatten([select_all(CompareSnvSens.val_by_freq_tsv),
                              select_all(CompareIndelSens.val_by_freq_tsv),
                              select_all(CompareSvSens.val_by_freq_tsv)]),
      n_key_columns = 3,
      out_prefix = sens_prefix + ".concordance_by_af",
      g2c_analysis_docker = g2c_analysis_docker
  }


  output {
    File? common_snv_ppv_bed = CollapseCommonSnvPpvSites.merged_file
    File? common_snv_sens_bed = CollapseCommonSnvSensSites.merged_file
    File? common_indel_ppv_bed = CollapseCommonIndelPpvSites.merged_file
    File? common_indel_sens_bed = CollapseCommonIndelSensSites.merged_file
    File? common_sv_ppv_bed = CollapseCommonSvPpvSites.merged_file
    File? common_sv_sens_bed = CollapseCommonSvSensSites.merged_file
    File ppv_by_size = SumPpvBySize.merged_distrib
    File sensitivity_by_size = SumSensBySize.merged_distrib
    File ppv_by_freq = SumPpvByFreq.merged_distrib
    File sensitivity_by_freq = SumSensByFreq.merged_distrib
  }
}


# Preps a sites.bed file for site comparison by generating:
# 1. A strict "query" set, which only includes variants with POS within an 
# eval interval, >50% coverage by all eval intervals, and size within [min_size, max_size]
# 2. A lenient "ref" set, which includes any variant overlapping any eval_interval
# and inclusive of all variants with size within [min_size / 3, 3 * max_size]
task PrepSites {
  input {
    Array[File] beds
    Array[File] bed_idxs
    File eval_interval_bed

    Int min_size = 0
    Int max_size = 1000000000
    Float lenient_size_scalar = 3.0
    Float strict_interval_coverage = 0.5
    
    String prefix

    Int disk_gb = 25
    
    String g2c_analysis_docker
  }

  parameter_meta {
    beds: {
      localization_optional: true
    }
  }

  Int loose_min_size = floor(min_size / lenient_size_scalar)
  Int loose_max_size = ceil(lenient_size_scalar * max_size)

  command <<<
    set -eu -o pipefail

    # Link all bed indexes to pwd
    while read bed_idx; do
      ln -s $bed_idx .
    done < ~{write_lines(bed_idxs)}

    # Save header line for first BED
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
    tabix --only-header -p bed -R ~{eval_interval_bed} ~{beds[0]} > header.bed

    # Loop over each bed and use tabix to remotely extract only the variants needed
    while read bed_uri; do
      export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`
      tabix -p bed -R ~{eval_interval_bed} $bed_uri
    done < ~{write_lines(beds)} \
    | awk -v FS="\t" -v min_size=~{loose_min_size} -v max_size=~{loose_max_size} \
      '{ if ($7>=min_size && $7<=max_size) print }' \
    | sort -Vk1,1 -k2,2n -k3,3n | uniq \
    | cat header.bed - | bgzip -c \
    > ~{prefix}.ref.bed.gz

    # Further filter the lenient "ref" set to produce a strict "query" set
    if [ $( zcat ~{prefix}.ref.bed.gz | fgrep -v "#" | wc -l ) -gt 0 ]; then
      /opt/pancan_germline_wgs/scripts/qc/vcf_qc/enforce_strict_intervals.py \
        -i ~{prefix}.ref.bed.gz \
        -t ~{eval_interval_bed} \
        -f ~{strict_interval_coverage} \
        -m ~{min_size} \
        -M ~{max_size} \
        -o ~{prefix}.query.bed.gz
    else
      cp ~{prefix}.ref.bed.gz ~{prefix}.query.bed.gz
    fi
  >>>

  output {
    File ref_bed = "~{prefix}.ref.bed.gz"
    File query_bed = "~{prefix}.query.bed.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


task CompareSites {
  input {
    File query_bed
    File ref_bed
    
    File genome_file
    String mode
    Float common_af
    Int overlap_pad = 1
    Boolean degenerate = false
    
    String prefix

    Int? disk_gb
    Float mem_gb = 3.75
    Int n_cpu = 2

    String g2c_analysis_docker
  }

  Int default_disk_gb = ceil(2.5 * size([query_bed, ref_bed], "GB")) + 10
  Int use_disk_gb = select_first([disk_gb, default_disk_gb])
  String degen_cmd = if degenerate then "--one-to-many" else ""

  command <<<
    set -eu -o pipefail

    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/compare_sites.py \
      -a ~{query_bed} \
      -b ~{ref_bed} \
      -g ~{genome_file} \
      -o ~{prefix} \
      --common-af ~{common_af} \
      --mode ~{mode} \
      --overlap-pad ~{overlap_pad} \
      ~{degen_cmd} \
      --no-reverse \
      --gzip
  >>>

  output {
    File all_sites_bed = "~{prefix}.loj.sites.bed.gz"
    File all_sites_bed_idx = "~{prefix}.loj.sites.bed.gz.tbi"
    File common_sites_bed = "~{prefix}.loj.sites.common.bed.gz"
    File common_sites_bed_idx = "~{prefix}.loj.sites.common.bed.gz.tbi"
    File val_by_size_tsv = "~{prefix}.loj.concordance_by_size.tsv.gz"
    File val_by_freq_tsv = "~{prefix}.loj.concordance_by_af.tsv.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk ~{use_disk_gb} HDD"
    preemptible: 3
  }
}
