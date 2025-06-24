# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Perform two-sided sample-level benchmarking of a "source" callset versus a reference "target" callset
# Note: this is a subroutine called within BenchmarkSamples.wdl, and is designed to be
# invoked by this parent workflow

# Expects that both callsets are provided according to G2C QC expectations


version 1.0


import "BenchmarkSites.wdl" as SiteBench
import "QcTasks.wdl" as QcTasks
import "Utilities.wdl" as Utils


workflow BenchmarkSamplesSingle {
  input {
    File source_all_sites_bed
    File? source_snv_bed
    File? source_indel_bed
    File? source_sv_bed
    File source_gt_tarball
    File source_sample_id_list
    String source_prefix

    File target_vcf
    File target_vcf_idx
    String target_prefix
    
    File id_map_tsv

    Array[File] eval_interval_beds
    Array[String]? eval_interval_bed_names
    File genome_file

    Int total_shards = 2500
    Int min_samples_per_shard = 10
    Float common_af_cutoff = 0.01

    String bcftools_docker
    String g2c_analysis_docker
  }

  # Subset target VCF to overlapping samples
  call QcTasks.MakeHeaderFiller {}
  call SubsetTargetVcf {
    input:
      vcf = target_vcf,
      vcf_idx = target_vcf_idx,
      supp_header = MakeHeaderFiller.supp_vcf_header,
      id_map_tsv = id_map_tsv,
      source_sample_id_list = source_sample_id_list,
      bcftools_docker = bcftools_docker
  }

  # Collect site metrics for target VCF
  call QcTasks.CollectSiteMetrics as TargetVcfToBed {
    input:
      vcf = SubsetTargetVcf.filtered_sites_vcf,
      vcf_idx = SubsetTargetVcf.filtered_sites_vcf_idx,
      n_samples = SubsetTargetVcf.n_samples_in_original_vcf,
      common_af_cutoff = common_af_cutoff,
      g2c_analysis_docker = g2c_analysis_docker
  }
  # TODO: concatenate target site metrics

  # Run site-level benchmarking to determine mapping of 
  # matching variant IDs between source & target
  call SiteBench.BenchmarkSites {
    input:
      source_snv_bed = source_snv_bed,
      source_indel_bed = source_indel_bed,
      source_sv_bed = source_sv_bed,
      source_prefix = source_prefix,
      target_snv_bed = TargetVcfToBed.snv_sites,
      target_indel_bed = TargetVcfToBed.indel_sites,
      target_sv_bed = TargetVcfToBed.sv_sites,
      target_prefix = target_prefix,
      eval_interval_beds = eval_interval_beds,
      eval_interval_bed_names = eval_interval_bed_names,
      genome_file = genome_file,
      total_shards = total_shards,
      common_af_cutoff = common_af_cutoff,
      make_full_id_maps = true,
      bcftools_docker = bcftools_docker,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Collect sample genotypes for subsetted target VCF
  call QcTasks.CollectSampleGenotypeMetrics as GatherTargetGts {
    input:
      vcf = SubsetTargetVcf.filtered_dense_vcf,
      vcf_idx = SubsetTargetVcf.filtered_dense_vcf_idx,
      site_metrics = TargetVcfToBed.all_sites,
      g2c_analysis_docker = g2c_analysis_docker
  }
  call QcTasks.ConcatGenotypeTsvs as MakeTargetGtTarball {
    input:
      tsvs = [GatherTargetGts.genotypes_tsv],
      output_prefix = target_prefix,
      g2c_analysis_docker = g2c_analysis_docker
  }

  # Shard the cleaned ID map to optimally parallelize GT benchmarking tasks
  Int n_tasks_per_shard = 2 * length(eval_interval_beds)
  Int n_naive_sample_splits = ceil(n_tasks_per_shard * SubsetTargetVcf.n_overlapping_samples / min_samples_per_shard)
  Int n_sample_shards = if n_naive_sample_splits < total_shards then n_naive_sample_splits else total_shards
  # call Utils.ShardTextFile as ShardSamples {
  #   input:
  #     input_file = SubsetTargetVcf.filtered_id_map,
  #     n_splits = n_sample_shards,
  #     out_prefix = "~{source_prefix}.~{target_prefix}.gt_bench.id_map_shard",
  #     shuffle = true,
  #     g2c_analysis_docker = g2c_analysis_docker
  # }

  # # Scatter over samples and perform benchmarking on each
  # scatter ( k in range(length(ShardSamples.shards)) ) {
  #   # TODO:
  #   # 5. Compare .tsvs from 4 to equivalent .tsvs generated from target/input callset while linking through variant ID map from #3
  #   call BenchmarkGenotypes as BenchGtPpv {
  #     input:
  #       variant_id_map = BenchmarkSites.ppv_variant_id_maps
  #   }

  #   # TODO: need to add sensitivity here
  # }

  # 6. Collapse outputs of 5 across all samples and compute summary metrics (medians?)

  output {}
}


# Subset a target VCF based on a dynamically-determined list of overlapping samples
task SubsetTargetVcf {
  input {
    File vcf
    File vcf_idx
    File supp_header
    File id_map_tsv
    File source_sample_id_list
    String bcftools_docker
  }

  Int disk_gb = ceil(2.5 * size(vcf, "GB")) + 10

  String filtered_map_fname = basename(id_map_tsv, ".tsv") + ".filtered.tsv"
  String filtered_dense_vcf_fname = basename(vcf, ".vcf.gz") + "filtered.vcf.gz"
  String filtered_sites_vcf_fname = basename(vcf, ".vcf.gz") + "filtered.sites.vcf.gz"

  command <<<
    set -eu -o pipefail

    # Get list of samples present in input VCF
    bcftools query -l ~{vcf} | sort > target.samples.list

    # Ensure all files are sorted for bash join
    sort ~{source_sample_id_list} > source.samples.list
    sort -k1,1 ~{id_map_tsv} > id_map.tsv

    # Use bash join to subset ID map
    join -j 1 -t $'\t' source.samples.list id_map.tsv \
    | sort -k2,2 \
    | join -1 2 -2 1 -t $'\t' - target.samples.list \
    | awk -v OFS="\t" '{ print $2, $1 }' \
    > ~{filtered_map_fname}

    # Subset input VCF to overlapping samples
    cut -f2 ~{filtered_map_fname} | sort -V | uniq > target.keep.list
    bcftools annotate -h ~{supp_header} ~{vcf} \
    | bcftools view --no-update --samples-file target.keep.list --force-samples \
    | bcftools view --no-update --include 'GT="alt" | FILTER="MULTIALLELIC"' \
      -Oz -o ~{filtered_dense_vcf_fname}
    tabix -p vcf -f ~{filtered_dense_vcf_fname}

    # Make a sites VCF for benchmarking comparison
    bcftools view --no-update -G \
      -Oz -o ~{filtered_sites_vcf_fname} \
      ~{filtered_dense_vcf_fname}
    tabix -p vcf -f ~{filtered_sites_vcf_fname}
  >>>

  output {
    File filtered_dense_vcf = filtered_dense_vcf_fname
    File filtered_dense_vcf_idx = "~{filtered_dense_vcf_fname}.tbi"

    File filtered_sites_vcf = filtered_sites_vcf_fname
    File filtered_sites_vcf_idx = "~{filtered_sites_vcf_fname}.tbi"

    File filtered_id_map = filtered_map_fname

    Int n_samples_in_original_vcf = length(read_lines("target.samples.list"))
    Int n_overlapping_samples = length(read_lines("target.keep.list"))
  }

  runtime {
    docker: bcftools_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


# Compare genotypes between callsets for a list of samples
task BenchmarkGenotypes {
  input {
    File variant_id_map
    File sample_id_map
    Boolean invert_sample_map = false
    String output_prefix

    File source_gt_tarball
    File target_gt_tarball

    File source_site_metrics
    Float common_af_cutoff = 0.01
    
    String g2c_analysis_docker
  }

  Int disk_gb = ceil(4 * size([source_gt_tarball, target_gt_tarball, source_site_metrics], "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Prep variant ID lists & metric files
    zcat ~{variant_id_map} | cut -f1 | sort -V | uniq > source.vids.list
    zcat ~{source_site_metrics} | head -n1 | fgrep "#" > site_metrics.header || true
    zcat ~{source_site_metrics} | fgrep -wf source.vids.list \
    | cat site_metrics.header - | bgzip -c \
    > source.metrics.bed.gz || true
    rm ~{source_site_metrics}
    zcat ~{variant_id_map} | cut -f2 | fgrep -xv "NA" | sort -V | uniq > target.vids.list

    # Prep sample lists
    if ~{invert_sample_map}; then
      awk -v FS="\t" -v OFS="\t" '{ print $2, $1 }' ~{sample_id_map} > sample.map.tsv
    else
      cp ~{sample_id_map} sample.map.tsv
    fi
    cut -f1 sample.map.tsv > source.samples.list
    cut -f2 sample.map.tsv > target.samples.list

    # Unpack source GTs and subset to samples of interest
    mkdir source_gts_raw
    tar -xzvf ~{source_gt_tarball} -C source_gts_raw/
    mkdir source_gts/
    while read sid; do
      find source_gts_raw/ -name "$sid.gt.sub.tsv.gz" \
      | xargs -I {} zcat {} | fgrep -wf source.vids.list \
      | gzip -c > source_gts/$sid.gt.tsv.gz || true
    done < source.samples.list
    rm -rf source_gts_raw ~{source_gt_tarball}

    # Unpack target GTs and subset to samples of interest
    mkdir target_gts_raw
    tar -xzvf ~{target_gt_tarball} -C target_gts_raw/
    mkdir target_gts/
    while read sid; do
      find target_gts_raw/ -name "$sid.gt.sub.tsv.gz" \
      | xargs -I {} zcat {} | fgrep -wf target.vids.list \
      | gzip -c > target_gts/$sid.gt.tsv.gz || true
    done < target.samples.list
    rm -rf target_gts_raw ~{target_gt_tarball}

    # Benchmark genotypes
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/compare_genotypes.R \
      --variant-map ~{variant_id_map} \
      --source-site-metrics source.metrics.bed.gz \
      --sample-map sample.map.tsv \
      --source-gt-dir source_gts/ \
      --target-gt-dir target_gts/ \
      --gt-tsv-suffix ".gt.sub.tsv.gz" \
      --common-af ~{common_af_cutoff} \
      --out-prefix ~{output_prefix}
    gzip -f ~{output_prefix}.gt_comparison.distrib.tsv
  >>>

  output {
    File gt_bench_distrib = "~{output_prefix}.gt_comparison.distrib.tsv.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}

