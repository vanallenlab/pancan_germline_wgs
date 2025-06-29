# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Generic WDL tasks used for VCF quality control


version 1.0


# Compare genotypes between callsets for a list of samples
task BenchmarkGenotypes {
  input {
    File variant_id_map
    File sample_id_map
    Boolean invert_sample_map = false
    File source_site_metrics
    String output_prefix

    File source_gt_tarball
    File target_gt_tarball

    Boolean report_by_gt = false
    Float common_af_cutoff = 0.01
    
    String g2c_analysis_docker
  }

  parameter_meta {
    source_gt_tarball: {
      localization_optional: true
    }
    target_gt_tarball: {
      localization_optional: true
    }
  }

  String gt_report_cmd = if report_by_gt then "--report-by-genotype" else ""
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
    zcat ~{variant_id_map} | cut -f2 | fgrep -xv "NA" | sort -V | uniq > target.vids.list || true

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
    gsutil -m cp ~{source_gt_tarball} source_gt_tarball.tar.gz
    tar -xzvf source_gt_tarball.tar.gz -C source_gts_raw/
    mkdir source_gts/
    while read sid; do
      find source_gts_raw/ -name "$sid.gt.tsv.gz" \
      | xargs -I {} zcat {} | fgrep -wf source.vids.list \
      | gzip -c > source_gts/$sid.gt.sub.tsv.gz || true
    done < source.samples.list
    rm -rf source_gts_raw source_gt_tarball.tar.gz
    echo "Contents of source_gts:"
    ls -lh source_gts/

    # Unpack target GTs and subset to samples of interest
    mkdir target_gts_raw
    gsutil -m cp ~{target_gt_tarball} target_gt_tarball.tar.gz
    tar -xzvf target_gt_tarball.tar.gz -C target_gts_raw/
    mkdir target_gts/
    while read sid; do
      find target_gts_raw/ -name "$sid.gt.tsv.gz" \
      | xargs -I {} zcat {} | fgrep -wf target.vids.list \
      | gzip -c > target_gts/$sid.gt.sub.tsv.gz || true
    done < target.samples.list
    rm -rf target_gts_raw target_gt_tarball.tar.gz
    echo "Contents of target_gts:"
    ls -lh source_gts/

    # Benchmark genotypes
    echo "Now benchmarking..."
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/compare_genotypes.R \
      --variant-map ~{variant_id_map} \
      --source-site-metrics source.metrics.bed.gz \
      --sample-map sample.map.tsv \
      --source-gt-dir source_gts/ \
      --target-gt-dir target_gts/ \
      --gt-tsv-suffix ".gt.sub.tsv.gz" \
      ~{gt_report_cmd} \
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


task CollectSampleGenotypeMetrics {
  input {
    File vcf
    File vcf_idx

    File? site_metrics
    Float? common_af_cutoff

    String g2c_analysis_docker
  }

  String out_base = basename(vcf, ".vcf.gz")
  String gt_outfile = out_base + ".genotypes.tsv.gz"
  String distrib_outfile = out_base + ".gt_distrib.tsv.gz"

  String site_metrics_basename = if defined(site_metrics) then basename(select_first([site_metrics])) else ""
  String distrib_cmd = if defined(site_metrics) then "--site-metrics ~{site_metrics_basename} --distrib-out ~{out_base}.gt_distrib.tsv.gz" else ""
  String common_cmd = if defined(site_metrics) && defined(common_af_cutoff) then "--common-af ~{common_af_cutoff}" else ""

  Int disk_gb = ceil(1.2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    if ~{defined(site_metrics)}; then
      mv ~{select_first([site_metrics])} ./
    fi

    bcftools query -i 'GT="alt"' -f '[%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVLEN\t%SAMPLE\t%GT\t%RD_CN\n]' ~{vcf} \
    | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/clean_sample_genotypes.py ~{distrib_cmd} ~{common_cmd} \
    | gzip -c \
    > ~{gt_outfile}
  >>>

  output {
    File genotypes_tsv = gt_outfile
    File? compressed_gt_distrib = distrib_outfile
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


task CollectSiteMetrics {
  input {
    File vcf
    File vcf_idx

    Int n_samples

    Float? min_af_bin
    Float? common_af_cutoff

    String g2c_analysis_docker
  }

  String out_prefix = basename(vcf, "vcf.gz")
  Int disk_gb = ceil(2 * size(vcf, "GB")) + 10

  String min_af_cmd = if defined(min_af_bin) then "--min-af-bin ~{min_af_bin}" else ""
  String common_cmd = if defined(common_af_cutoff) then "--common-af ~{common_af_cutoff}" else ""

  command <<<
    set -eu -o pipefail

    # Ensure all necessary fields are defined in VCF header
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">" > header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_FREQ,Number=1,Type=Float,Description=\"CNV frequency\">" >> header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_COUNT,Number=1,Type=Integer,Description=\"Nondip count.\">" >> header.supp.vcf
    echo "##INFO=<ID=HWE,Number=A,Type=Float,Description=\"HWE test\">" >> header.supp.vcf
    echo "##INFO=<ID=ExcHet,Number=A,Type=Float,Description=\"ExcHet test\">" >> header.supp.vcf

    # Collect stats and split into SNV, indel, and SV files
    bcftools annotate -h header.supp.vcf ~{vcf} \
    | bcftools query \
      -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%INFO/SVLEN\t%INFO/AN\t%INFO/AC\t%INFO/AF\t%INFO/CN_NONREF_COUNT\t%INFO/CN_NONREF_FREQ\t%INFO/AC_Het\t%INFO/AC_Hom\t%INFO/AC_Hemi\t%INFO/HWE\n' \
    | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/clean_site_metrics.py \
      ~{min_af_cmd} \
      ~{common_cmd} \
      -o ~{out_prefix} \
      --gzip \
      -N ~{n_samples}

    # Concatenate all site metrics for downstream 
    # compatability with CollectSampleGenotypeMetrics
    find ./ -name "~{out_prefix}.*.sites.bed.gz" > site_beds.list
    zcat $( head -n1 site_beds.list ) | head -n1 > site_metrics.header || true
    while read site_file; do
      zcat $site_file | fgrep -v "#" || true
    done < site_beds.list \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
    | uniq \
    | cat site_metrics.header - \
    | bgzip -c \
    > ~{out_prefix}.all.sites.bed.gz || true
    tabix -p bed -f ~{out_prefix}.all.sites.bed.gz
  >>>

  output {
    File all_sites = out_prefix + ".all.sites.bed.gz"
    File all_sites_idx = out_prefix + ".all.sites.bed.gz.tbi"
    File? snv_sites = out_prefix + ".snv.sites.bed.gz"
    File? snv_sites_idx = out_prefix + ".snv.sites.bed.gz.tbi"
    File? indel_sites = out_prefix + ".indel.sites.bed.gz"
    File? indel_sites_idx = out_prefix + ".indel.sites.bed.gz.tbi"
    File? sv_sites = out_prefix + ".sv.sites.bed.gz"
    File? sv_sites_idx = out_prefix + ".sv.sites.bed.gz.tbi"
    File? common_snv_sites = out_prefix + ".snv.sites.common.bed.gz"
    File? common_snv_sites_idx = out_prefix + ".snv.sites.common.bed.gz.tbi"
    File? common_indel_sites = out_prefix + ".indel.sites.common.bed.gz"
    File? common_indel_sites_idx = out_prefix + ".indel.sites.common.bed.gz.tbi"
    File? common_sv_sites = out_prefix + ".sv.sites.common.bed.gz"
    File? common_sv_sites_idx = out_prefix + ".sv.sites.common.bed.gz.tbi"
    File size_distrib = out_prefix + ".size_distrib.tsv.gz"
    File af_distrib = out_prefix + ".af_distrib.tsv.gz"
    File size_vs_af_distrib = out_prefix + ".size_vs_af_distrib.tsv.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


# Concatenate an array of sample genotype .tsvs and split into one file per sample
task ConcatGenotypeTsvs {
  input {
    Array[File] tsvs
    String output_prefix
    String g2c_analysis_docker
  }

  Int disk_gb = ceil(2 * size(tsvs, "GB")) + 10
  String outdir = output_prefix + "_sample_genotypes"

  command <<<
    set -eu -o pipefail

    mkdir ~{outdir}

    # Split each input tsv by appending to sample-specific files in outdir/
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/split_and_append_gts.py \
      ~{write_lines(tsvs)} \
      ~{outdir}

    # Compress outdir as tarball for easier GCP/Cromwell IO
    tar -czvf ~{outdir}.tar.gz ~{outdir}
  >>>

  output {
    File genotypes_tarball = "~{outdir}.tar.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}


# Make an empty text file as a placeholder
task MakeDummyFile {
  input {}

  command <<<
    touch "dummy.txt"
  >>>

  output {
    File empty_file = "dummy.txt"
  }

  runtime {
    docker: "marketplace.gcr.io/google/ubuntu1804"
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 15 HDD"
    preemptible: 3
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


# Divide a BED-style intervals file into shards balanced on total genomic footprint
task ShardIntervals {
  input {
    File intervals_bed
    Int n_shards
    String prefix
    String g2c_analysis_docker
  }

  Int default_disk_gb = ceil(10 * size(intervals_bed, "GB")) + 15

  command <<<
    set -eu -o pipefail

    # Ensure intervals file is simple, sorted, bgzipped, and merged
    zcat ~{intervals_bed} \
    | fgrep -v "#" \
    | cut -f1-3 \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    | bgzip -c \
    > intervals.prepped.bed.gz

    # Estimate total size per shard
    shard_size=$( zcat intervals.prepped.bed.gz \
                  | awk -v denom=~{n_shards} \
                    '{ sum+=$3-$2 }END{ printf "%.0f", sum / denom }' )

    # Ensure no interval is larger than shard_size
    /opt/pancan_germline_wgs/scripts/gatkhc_helpers/split_intervals.py \
      -i intervals.prepped.bed.gz \
      -t $shard_size \
      --max-interval-size $shard_size \
      --bed-style \
      -o intervals.prepped.split.bed
    bgzip intervals.prepped.split.bed

    # Shard intervals
    k=1
    bp_remain=$shard_size
    while read chrom start end; do

      size=$(( $end - $start ))

      if [ $size -le $bp_remain ]; then
        echo -e "$chrom\t$start\t$end" >> "~{prefix}.$k.bed"
        bp_remain=$(( $bp_remain - $size ))

      else
        new_end=$(( $start + $bp_remain ))
        echo -e "$chrom\t$start\t$new_end" >> "~{prefix}.$k.bed"
        
        # Because we split all intervals to be <= shard_size above,
        # we know that the remainder must always be < shard_size
        ((k++))
        start=$new_end
        size=$(( $end - $start ))
        echo -e "$chrom\t$start\t$end" >> "~{prefix}.$k.bed"
        bp_remain=$(( $shard_size - $size ))
      fi

    done < <( zcat intervals.prepped.split.bed.gz )
    
    find ./ -name "~{prefix}.*.bed" | xargs -I {} bgzip {}
  >>>

  output {
    Array[File] interval_shards = glob("~{prefix}.*.bed.gz")
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk ~{default_disk_gb} HDD"
    preemptible: 3
  }
}


# Make an empty benchmark bed to backfill for missing values in optional arrays
task MakeEmptyBenchBed {
  input {
    String docker
    String header_format = "site"
  }

  command <<<
    set -eu -o pipefail

    if [ "~{header_format}" == "site" ]; then
      echo -e "#chrom\tstart\tend\tvid\tclass\tsubclass\tsize\tac\taf\tfreq_het\tfreq_hom\thwe" \
      | bgzip -c > dummy.bed.gz
    fi
  >>>

  output {
    File empty_bed = "dummy.bed.gz"
  }

  runtime {
    docker: docker
    memory: "1.5 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
  }
}


task MakeHeaderFiller {
  input {}

  command <<<
    set -eu -o pipefail

    echo "##FILTER=<ID=MULTIALLELIC,Description=\"Multiallelic site\">" > header.supp.vcf
    echo "##FORMAT=<ID=RD_CN,Number=1,Type=Integer,Description=\"Predicted copy state\">" >> header.supp.vcf
    echo "##INFO=<ID=AC_Het,Number=A,Type=Integer,Description=\"Heterozygous allele counts\">" >> header.supp.vcf
    echo "##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description=\"Homozygous allele counts\">" >> header.supp.vcf
    echo "##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description=\"Hemizygous allele counts\">" >> header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_COUNT,Number=1,Type=Integer,Description=\"Nondip count.\">" >> header.supp.vcf
    echo "##INFO=<ID=CN_NONREF_FREQ,Number=1,Type=Float,Description=\"CNV frequency\">" >> header.supp.vcf
    echo "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">" >> header.supp.vcf
    echo "##INFO=<ID=HWE,Number=A,Type=Float,Description=\"HWE test\">" >> header.supp.vcf
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">" >> header.supp.vcf
    echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">" >> header.supp.vcf
    
  >>>

  output {
    File supp_vcf_header = "header.supp.vcf"
  }

  runtime {
    docker: "marketplace.gcr.io/google/ubuntu1804"
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
    max_retries: 1
  }
}


# Sum one or more compressed distributions generated by CollectSiteMetrics
task SumCompressedDistribs {
  input {
    Array[File] distrib_tsvs
    Int n_key_columns = 2
    String out_prefix
    String g2c_analysis_docker
  }

  Int disk_gb = ceil(1.5 * size(distrib_tsvs, "GB")) + 10

  command <<<
    set -eu -o pipefail
    
    /opt/pancan_germline_wgs/scripts/qc/vcf_qc/sum_compressed_distribs.py \
      -o "~{out_prefix}.merged.tsv" \
      -k ~{n_key_columns} \
      ~{sep=" " distrib_tsvs}
    gzip -f "~{out_prefix}.merged.tsv"
  >>>

  output {
    File merged_distrib = "~{out_prefix}.merged.tsv.gz"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "1.75 GB"
    cpu: 1
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
  }
}

