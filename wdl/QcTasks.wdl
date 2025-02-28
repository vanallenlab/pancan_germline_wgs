# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Generic WDL tasks used for VCF quality control


version 1.0


task CollectSiteMetrics {
  input {
    File vcf
    File vcf_idx

    Int n_samples

    String g2c_analysis_docker
  }

  String out_prefix = basename(vcf, ".vcf.gz")
  Int disk_gb = ceil(2 * size(vcf, "GB")) + 10

  command <<<
    set -eu -o pipefail

    # Ensure all necessary fields are defined in VCF header
    echo "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length\">" > header.supp.vcf
    echo "##INFO=<ID=HWE,Number=A,Type=Float,Description=\"HWE test\">" >> header.supp.vcf
    echo "##INFO=<ID=ExcHet,Number=A,Type=Float,Description=\"ExcHet test\">" >> header.supp.vcf

    # Collect stats and split into SNV, indel, and SV files
    bcftools annotate -h header.supp.vcf ~{vcf} \
    | bcftools query \
      -f '%CHROM\t%POS\t%END\t%REF\t%ALT\t%INFO/SVLEN\t%INFO/AC\t%INFO/AF\t%INFO/HWE\t%INFO/ExcHet\n' \
    | /opt/pancan_germline_wgs/scripts/qc/vcf_qc/clean_site_metrics.py \
      -o ~{out_prefix} \
      --gzip \
      -N ~{n_samples}

  >>>

  output {
    # TBD
  }

  runtime {
    docker: g2c_analysis_docker
    memory: "3.5 GB"
    cpu: 2
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    max_retries: 1
  }
}
