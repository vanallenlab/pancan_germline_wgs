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
    File genome_file

    Float recip_overlap = 0.98
    Int bp_dist = 500
    Float sample_overlap = 0.05
    String cluster_sv_types = "DEL,DUP,INS,INV,CPX"
    Float min_af = 0.001
    Int min_ac = 10

    String g2c_analysis_docker
  }

  call DefineClusters {
    input:
      vcf = vcf,
      vcf_idx = vcf_idx,
      genome_file = genome_file,
      recip_overlap = recip_overlap,
      bp_dist = bp_dist,
      sample_overlap = sample_overlap,
      cluster_sv_types = cluster_sv_types,
      min_af = min_af,
      min_ac = min_ac,
      g2c_analysis_docker = g2c_analysis_docker    
  }

  if ( DefineClusters.n_clusters > 0 ) {
    call ResolveClusters {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        clusters = DefineClusters.clusters,
        g2c_analysis_docker = g2c_analysis_docker
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
    File genome_file

    Float recip_overlap
    Int bp_dist
    Float sample_overlap
    String cluster_sv_types
    Float min_af
    Int min_ac
    
    Float mem_gb = 7.5
    Int n_cpu = 4

    String g2c_analysis_docker
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
    echo "input.vcf.gz" > input_vcf.list
    svtk vcfcluster \
      -d ~{bp_dist} \
      -f ~{recip_overlap} \
      -o ~{sample_overlap} \
      -p ~{out_prefix}.reclustered \
      -t "~{cluster_sv_types}" \
      --preserve-header \
      --skip-merge \
      input_vcf.list \
      - \
    | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%INFO/CLUSTER\n' \
    > ~{out_prefix}.sv_cluster_assignments.bed

    # Second pass for clusters of two or more strictly identical variants
    # This step is not filtered by frequency, and is pre-screened for positions
    # where at least two variants share the same start position
    bcftools query -f '%CHROM\t%POS\n' ~{vcf} \
    | sort -Vk1,1 -k2,2n | uniq -c \
    | awk -v OFS="\t" '{ if ($1>1) print $2, $3-2, $3+2 }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    > candidate_identical_loci.bed
    /opt/pancan_germline_wgs/scripts/utilities/filter_vcf_by_pos.py \
      -i ~{vcf} \
      -r candidate_identical_loci.bed \
      -g ~{genome_file} \
    | bcftools annotate \
      -h add_header.vcf \
      -Oz -o input2.vcf.gz
    tabix -f input2.vcf.gz
    echo "input2.vcf.gz" > input2_vcf.list
    svtk vcfcluster \
      -d 1 \
      -f 1 \
      -p "~{out_prefix}.identical" \
      -t "~{cluster_sv_types}" \
      --preserve-header \
      --skip-merge \
      input2_vcf.list \
      - \
    | bcftools query -f '%CHROM\t%POS\t%END\t%ID\t%INFO/CLUSTER\n' \
    > ~{out_prefix}.identical_assignment.bed
    cut -f5 ~{out_prefix}.identical_assignment.bed \
    | uniq -c | awk '{ if ($1>1) print $2 }' \
    > candidate_cluster_idxs.list
    while read cidx; do
      awk -v cidx=$cidx -v OFS="\t" \
        '{ if ($5==cidx) print }' \
        ~{out_prefix}.identical_assignment.bed
    done < candidate_cluster_idxs.list \
    > ~{out_prefix}.identical_assignment.clusters.bed

    # Update rough cluster assignments to ensure that 
    # all strictly identical variants are linked. Outputs
    # a BED4 file with one row per multi-variant cluster
    /opt/pancan_germline_wgs/scripts/gatksv_helpers/update_cluster_assignments.py \
      -i ~{out_prefix}.sv_cluster_assignments.bed \
      -u ~{out_prefix}.identical_assignment.clusters.bed \
      -o ~{out_prefix}.sv_clusters.bed
  >>>

  output {
    File clusters = "~{out_prefix}.sv_clusters.bed"
    Int n_clusters = length(read_lines("~{out_prefix}.sv_clusters.bed"))
  }

  runtime {
    docker: g2c_analysis_docker
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

    String g2c_analysis_docker
  }

  String out_prefix = basename(vcf, ".vcf.gz")
  String vcf_out =  "~{out_prefix}.reclustered.vcf.gz"

  Int disk_gb = ceil(10 * size(vcf, "GB")) + 10

  Float sort_mem = 2 * mem_gb / 3

  command <<<
    set -eu -o pipefail

    # Process each variant cluster and write to reclustered output file
    /opt/pancan_germline_wgs/scripts/gatksv_helpers/resolve_variant_clusters.py \
      --in-vcf ~{vcf} \
      --clusters ~{clusters} \
      --prefix "~{out_prefix}_reclustered" \
    | bcftools sort -m "~{sort_mem}G" \
    | bcftools +fill-tags -Oz -o ~{out_prefix}.resolved_clusters.vcf.gz -- -t AC,AN,AF
    tabix -p vcf -f ~{out_prefix}.resolved_clusters.vcf.gz

    # Remove cluster members from input VCF and merge remainder with new clustered records
    cut -f4 ~{clusters} | sed 's/,/\n/g' | sort | uniq > member.vids.list
    bcftools view --exclude 'ID=@member.vids.list' ~{vcf} -Oz -o ~{out_prefix}.remainder.vcf.gz
    tabix -p vcf -f ~{out_prefix}.remainder.vcf.gz
    bcftools concat \
      --allow-overlaps \
      ~{out_prefix}.remainder.vcf.gz \
      ~{out_prefix}.resolved_clusters.vcf.gz \
      -Oz -o ~{vcf_out} 
    tabix -p vcf -f ~{vcf_out}
  >>>

  output {
    File reclustered_vcf = vcf_out
    File reclustered_vcf_idx = "~{vcf_out}.tbi"
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "~{mem_gb} GB"
    cpu: n_cpu
    disks: "local-disk ~{disk_gb} HDD"
    preemptible: 3
    maxRetries: 1
    bootDiskSizeGb: 30
  }
}
