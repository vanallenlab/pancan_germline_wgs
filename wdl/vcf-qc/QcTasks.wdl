# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Generic WDL tasks used for VCF quality control


version 1.0


task CollectSampleGenotypeMetrics {
  input {
    File vcf
    File vcf_idx
    File? site_metrics

    String g2c_analysis_docker
  }

  String out_base = basename(vcf, ".vcf.gz")
  String gt_outfile = out_base + ".genotypes.tsv.gz"
  String distrib_outfile = out_base + ".gt_distrib.tsv.gz"

  String distrib_cmd = if defined(site_metrics) then "--site-metrics ~{select_first([site_metrics])} --distrib-out ~{out_base}.gt_distrib.tsv.gz" else ""

  Int disk_gb = ceil(1.2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    bcftools query -i 'GT="alt"' -f '[%SAMPLE\t%ID\t%GT\t%RD_CN\n]' ~{vcf} \
    | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/clean_sample_genotypes.py ~{distrib_cmd} \
    | gzip -c \
    > ~{gt_outfile}
  >>>

  output {
    File genotypes_tsv = gt_outfile
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

  String out_prefix = basename(vcf, ".vcf.gz")
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
    find ~{out_prefix} -name "*.sites.bed.gz" > site_beds.list
    head -n1 $( head -n1 site_beds.list ) > site_metrics.header
    while read site_file; do
      zcat $site_file | fgrep -v "#"
    done < site_beds.list \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
    | uniq \
    | cat site_metrics.header \
    | bgzip -c \
    > ~{out_prefix}.all.sites.bed.gz
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

