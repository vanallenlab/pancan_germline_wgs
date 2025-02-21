# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Post hoc normalization & variant counting from a GATK-HC joint genotyped VCF

# Expected input follows the conventions of GnarlyJointGenotypingPart1.wdl


version 1.0


import "Utilities.wdl" as Utils


workflow PosthocCleanupPart1 {
  input {
    File vcf
    File vcf_idx

    File unpadded_intervals_file
    File ref_fasta
    
    String linux_docker
    String bcftools_docker
    String g2c_pipeline_docker
  }

  # Parse intervals file for scattering
  call Utils.ParseIntervals {
    input:
      intervals_list = unpadded_intervals_file,
      docker = linux_docker
  }

  # Parallelize over scatter intervals
  scatter ( interval_info in ParseIntervals.interval_info ) {

    String interval_coords = interval_info.right
    String shard_prefix = basename(vcf, ".vcf.gz") + "." + interval_info.left

    # Additional cleanup step added for G2C to coerce to parsimonious format
    call NormalizeShortVariants as NormalizeVcf {
      input:
        vcf = vcf,
        vcf_idx = vcf_idx,
        ref_fasta = ref_fasta,
        interval = interval_coords,
        bcftools_docker = bcftools_docker
    }

    # Count number of non-reference SNVs, insertions, and deletions per sample
    call CountShortVariantsPerSample {
      input:
        vcf = NormalizeVcf.norm_vcf,
        vcf_idx = NormalizeVcf.norm_vcf_idx,
        bcftools_docker = bcftools_docker
    }
  }

  # Combine sharded VCFs
  call Utils.ConcatVcfs {
    input:
      vcfs = NormalizeVcf.norm_vcf,
      vcf_idxs = NormalizeVcf.norm_vcf_idx,
      out_prefix = basename(vcf, "vcf.gz") + ".norm.vcf.gz"
  }

  # Sum variant counts
  call Utils.SumSvCountsPerSample as SumCounts {
    input:
      count_tsvs = CountShortVariantsPerSample.counts_tsv,
      output_prefix = basename(vcf, ".vcf.gz") + ".norm",
      docker = g2c_pipeline_docker
  }

  output {
    File normalized_vcf = ConcatVcfs.merged_vcf
    File normalized_vcf_idx = ConcatVcfs.merged_vcf_idx
    File counts_per_sample = SumCounts.summed_tsv
  }
}


task CountShortVariantsPerSample {
  input {
    File vcf
    File vcf_idx
    String bcftools_docker
  }

  String outfile = basename(vcf, ".vcf.gz") + ".variant_counts.tsv"
  Int disk_gb = ceil(1.2 * size([vcf], "GB")) + 10

  command <<<
    set -eu -o pipefail

    echo -e "sample\ttype\tcount" > ~{outfile}

    bcftools query \
      --include 'GT == "alt"' \
      --format '[%SAMPLE\t%REF\t%ALT\n]' \
      ~{vcf} \
    | awk -v FS="\t" -v OFS="\t" \
      '{ vl=length($3)-length($2); \
         if (vl>0) { vt="INS" }else \
         if (vl<0) { vt="DEL" }else
         { vt="SNV" };
         print $1, vt }' \
    | sort -Vk1,1 -k2,2V | uniq -c \
    | awk -v OFS="\t" '{ print $2, $3, $1 }' \
    >> ~{outfile}
  >>>

  output {
    File counts_tsv = outfile
  }

  runtime {
    docker: bcftools_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}


task NormalizeShortVariants {
  input {
    File vcf
    File vcf_idx
    File ref_fasta
    String interval
    String bcftools_docker
  }

  String outfile = basename(vcf, ".vcf.gz") + ".norm.vcf.gz"
  Int disk_gb = ceil(3 * size([vcf, ref_fasta], "GB")) + 10

  parameter_meta {
    vcf: {
      localization_optional: true
    }
  }

  command <<<
    set -eu -o pipefail

    # Symlink vcf_idx to current working dir
    ln -s ~{vcf_idx} .

    # Enable remote streaming of VCF
    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    bcftools norm \
      --regions "~{interval}" \
      --fasta-ref ~{ref_fasta} \
      --check-ref s \
      --multiallelics - \
      --site-win 10000 \
      --threads 4 \
      ~{vcf} \
    | bcftools view \
      --exclude 'alt[0] == "*"' \
      -Oz -o ~{outfile}

    tabix -p vcf ~{outfile}
  >>>

  output {
    File norm_vcf = outfile
    File norm_vcf_idx = "~{outfile}.tbi"
  }

  runtime {
    docker: bcftools_docker
    memory: "7.5 GB"
    cpu: 4
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    maxRetries: 1
  }
}