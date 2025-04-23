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


workflow BenchmarkSitesSingle {
  input {
    File eval_interval_bed
    String eval_prefix

    Pair[File, File]? source_snv_bed_and_idx
    Pair[File, File]? target_snv_bed_and_idx

    Pair[File, File]? source_indel_bed_and_idx
    Pair[File, File]? target_indel_bed_and_idx

    Pair[File, File]? source_sv_bed_and_idx
    Pair[File, File]? target_sv_bed_and_idx

    String source_prefix
    String target_prefix

    Int total_shards

    String g2c_analysis_docker
  }

  Boolean do_snv = if defined(source_snv_bed_and_idx) 
                   && defined(target_snv_bed_and_idx)
                   then true else false
  Boolean do_indel = if defined(source_indel_bed_and_idx) 
                     && defined(target_indel_bed_and_idx) 
                     then true else false
  Boolean do_sv = if defined(source_sv_bed_and_idx) 
                  && defined(target_sv_bed_and_idx) 
                  then true else false
  Int n_var_types = (if do_snv then 1 else 0) + 
                    (if do_indel then 1 else 0) + 
                    (if do_sv then 1 else 0)

  # Shard evaluation intervals, balancing on genomic bp
  Int target_shards = floor(total_shards / n_var_types)
  call QcTasks.ShardIntervals {
    input:
      intervals_bed = eval_interval_bed,
      n_shards = target_shards,
      prefix = eval_prefix,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Process each eval interval in parallel
  scatter ( eval_shard in ShardIntervals.interval_shards ) {
    # Compare SNVs, if provided
    # TODO: implement this

    # Compare indels, if provided
    # TODO: implement this
    # Note that we need to consider small SVs here too

    # Compare SVs, if provided
    # Note that we need to consider large indels here too
  }

  # Collapse SNV comparison results
  # TODO: implement this

  # Collapse indel comparison results
  # TODO: implement this

  # Collapse SV comparison results
  # TODO: implement this

  output {}
}

