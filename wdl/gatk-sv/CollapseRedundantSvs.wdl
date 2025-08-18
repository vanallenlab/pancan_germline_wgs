# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Collapse redundant SVs at the end of the GATK-SV pipeline
# Designed to be implemented after module 19 (FilterGenotypes)


version 1.0


workflow CollapseRedundantSvs {
  input {
    File vcf
    File vcf_idx

    Float recip_overlap = 0.98
    Int bp_dist = 500
    Float sample_overlap = 0.05
    String cluster_sv_types = "DEL,DUP,INS,INV,CPX"
    Float min_af = 0.001
    Int min_ac = 10

    String g2c_pipeline_docker
  }

  call DefineClusters {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      recip_overlap = recip_overlap,
      bp_dist = bp_dist,
      sample_overlap = sample_overlap,
      cluster_sv_types = cluster_sv_types,
      min_af = min_af,
      min_ac = min_ac,
      g2c_pipeline_docker = g2c_pipeline_docker    
  }

  if ( DefineClusters.n_clusters > 0 ) {
    call ResolveClusters {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        clusters = DefineClusters.clusters,
        g2c_pipeline_docker = g2c_pipeline_docker
    }
  }

  output {
    File reclustered_vcf = select_first([ResolveClusters.reclustered_vcf, vcf])
    File reclustered_vcf_idx = select_first([ResolveClusters.reclustered_vcf_idx, vcf_idx])
  }
}


task DefineClusters {
  input {
    File vcf
    File vcf_idx

    Float recip_overlap
    Int bp_dist
    Float sample_overlap
    String cluster_sv_types
    Float min_af
    Int min_ac
    
    Float mem_gb = 7.5
    Int n_cpu = 4

    String g2c_pipeline_docker
  }

  String out_prefix = basename(vcf, ".vcf.gz")

  Int disk_gb = ceil(2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Assign variants to clusters using svtk
    echo "##source=gatksv" > add_header.vcf
    bcftools annotate -h add_header.vcf ~{vcf} \
    | bcftools +fill-tags -- -t AF,AC \
    | bcftools view --include "AF>=~{min_af} & AC>=~{min_ac}" -Oz -o input.vcf.gz
    tabix -f input.vcf.gz
    svtk vcfcluster \
      -d ~{bp_dist} \
      -f ~{recip_overlap} \
      -o ~{sample_overlap} \
      -p "~{out_prefix}.reclustered" \
      -t "~{cluster_sv_types}" \
      --preserve-header \
      --skip-merge \
      <( echo "input.vcf.gz" ) \
      - \
    | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%INFO/CLUSTER\n' \
    > ~{out_prefix}.sv_cluster_assignments.tsv

    # Identify clusters of two or more variants
    while read cidx; do
      awk -v cidx=$cidx -v OFS="\t" \
        '{ if ($5==cidx) print $1, $2, $3, $4 }' \
        ~{out_prefix}.sv_cluster_assignments.tsv \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | bedtools merge -i - -d ~{bp_dist} -c 4 -o distinct
    done < <( cut -f5 ~{out_prefix}.sv_cluster_assignments.tsv \
              | uniq -c | awk '{ if ($1>1) print $2 }' ) \
    | awk '{ if ($4~/,/) print }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    > ~{out_prefix}.sv_clusters.bed
  >>>

  output {
    File clusters = "~{out_prefix}.sv_clusters.bed"
    Int n_clusters = length(read_lines("~{out_prefix}.sv_clusters.bed"))
  }

  runtime {
    docker: g2c_pipeline_docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
    maxRetries: 1
  }
}



task ResolveClusters {
  input {
    File vcf
    File vcf_idx
    File clusters

    Float mem_gb = 7.5
    Int n_cpu = 4

    String g2c_pipeline_docker
  }

  String out_prefix = basename(vcf, ".vcf.gz")
  String vcf_out =  "~{out_prefix}.reclustered.vcf.gz"

  Int disk_gb = ceil(5 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Process each variant cluster and write to reclustered output file
    /opt/pancan_germline_wgs/scripts/gatksv_helpers/resolve_variant_clusters.py \
      --in-vcf ~{vcf} \
      --clusters ~{clusters} \
      --prefix ~{out_prefix} \
    | bcftools +fill-tags -Oz -o ~{out_prefix}.resolved_clusters.vcf.gz -- -t AC,AN,AF
    tabix -p vcf -f ~{out_prefix}.resolved_clusters.vcf.gz

    # Remove cluster members from input VCF and merge remainder with new clustered records
    cut -f4 ~{clusters} | sed 's/,/\n/g' | sort | uniq > member.vids.list
    bcftools view --exclude 'ID=@member.vids.list' ~{vcf} -Oz -o ~{out_prefix}.remainder.vcf.gz
    tabix -p vcf -f ~{out_prefix}.remainder.vcf.gz
    bcftools concat \
      --allow-overlaps \
      -Oz -o ~{vcf_out} \
      ~{out_prefix}.remainder.vcf.gz \
      ~{out_prefix}.resolved_clusters.vcf.gz
    tabix -p vcf -f ~{vcf_out}
  >>>

  output {
    File reclustered_vcf = vcf_out
    File reclustered_vcf_idx = "~{vcf_out}.tbi"
  }

  runtime {
    docker: g2c_pipeline_docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
    maxRetries: 1
  }
}
