# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform two-sided site-level benchmarking of a "source" callset versus a reference "target" callset

# Expects that both callsets have already been processed by CollectVcfQcMetrics.wdl


version 1.0


import "Utilities.wdl" as Utils
import "BenchmarkSitesSingle.wdl" as BenchSingle


workflow BenchmarkSites {
  input {
    File? source_snv_bed
    File? source_snv_bed_idx
    File? source_indel_bed
    File? source_indel_bed_idx
    File? source_sv_bed
    File? source_sv_bed_idx
    String source_prefix

    File? target_snv_bed
    File? target_snv_bed_idx
    File? target_indel_bed
    File? target_indel_bed_idx
    File? target_sv_bed
    File? target_sv_bed_idx
    String target_prefix

    Array[File] eval_interval_beds
    Array[String]? eval_interval_bed_names

    Int total_shards = 2500

    String bcftools_docker
    String g2c_analysis_docker
  }

  Boolean do_snv = if defined(source_snv_bed) && defined(target_snv_bed) then true else false
  Boolean do_indel = if defined(source_indel_bed) && defined(target_indel_bed) then true else false
  Boolean do_sv = if defined(source_sv_bed) && defined(target_sv_bed) then true else false

  # Index SNV BEDs as necessary
  if ( do_snv ){
    if ( !defined(source_snv_bed_idx) ){
      call Utils.IndexBed as IndexSourceSnvs {
        input:
          bed = select_first([source_snv_bed]),
          docker = bcftools_docker
      }
    }
    if ( !defined(target_snv_bed_idx) ){
      call Utils.IndexBed as IndexTargetSnvs {
        input:
          bed = select_first([target_snv_bed]),
          docker = bcftools_docker
      }
    }
  }

  # Index indel BEDs as necessary
  if ( do_indel ){
    if ( !defined(source_indel_bed_idx) ){
      call Utils.IndexBed as IndexSourceIndels {
        input:
          bed = select_first([source_indel_bed]),
          docker = bcftools_docker
      }
    }
    if ( !defined(target_indel_bed_idx) ){
      call Utils.IndexBed as IndexTargetIndels {
        input:
          bed = select_first([target_indel_bed]),
          docker = bcftools_docker
      }
    }
  }

  # Index SV BEDs as necessary
  if ( do_sv ){
    if ( !defined(source_sv_bed_idx) ){
      call Utils.IndexBed as IndexSourceSvs {
        input:
          bed = select_first([source_sv_bed]),
          docker = bcftools_docker
      }
    }
    if ( !defined(target_sv_bed_idx) ){
      call Utils.IndexBed as IndexTargetSvs {
        input:
          bed = select_first([target_sv_bed]),
          docker = bcftools_docker
      }
    }
  }

  Int n_eval_beds = length(eval_interval_beds)
  Int shards_per_eval_bed = floor(total_shards / n_eval_beds)

  scatter ( eval_idx in range(n_eval_beds) ){

    File eval_bed = eval_interval_beds[eval_idx]
    if ( defined(eval_interval_bed_names) ) {
      String eval_name_custom = select_first([eval_interval_bed_names])[eval_idx]
    }
    String eval_name = select_first([eval_name_custom, basename(eval_bed, ".bed.gz")])

    # Run evaluation for a single eval BED as a subworkflow
    call BenchSingle.BenchmarkSitesSingle as BenchmarkTask {

    }
    # TODO: implement this

  }


  output {}
}

