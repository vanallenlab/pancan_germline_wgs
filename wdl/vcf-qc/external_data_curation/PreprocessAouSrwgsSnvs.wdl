# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Preprocess existing AoU short-read WGS-based SNV/indel calls to match G2C QC expectations
# Note that this is expected to operate on a single chromosome (unlike the other prepreocessing WDLs)


version 1.0


import "QcTasks.wdl" as QcTasks


workflow PreprocessAouSrwgsSnvs {
  input {
	  String vds_uri = "gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/vds/hail.vds"
    File? samples_list

    String g2c_pipeline_docker
  }

  call QcTasks.MakeHeaderFiller {}

  output {}
}


task ExportIntervalVcf {
  input {
    String vds_uri
    File samples_list
    String interval
  }

  command <<<
    set -eu -o pipefail

    python3 << CODE
import hail as hl

hl.init(default_reference="GRCh38")

# Inputs
vds_uri = "${vds_uri}"
interval_str = "${interval}"
output_path = f"{output_bucket}/region_{interval_str.replace(':', '_').replace('-', '_')}.vcf.bgz"

# Parse interval string into a Hail interval
contig, positions = interval_str.split(":")
start, end = map(int, positions.split("-"))
interval = hl.parse_locus_interval(f"{contig}:{start}-{end}")

# Read and filter VDS
vds = hl.vds.read_vds(vds_uri)

# Sample list
samples_table = hl.import_table("${samples_list}", no_header=True).key_by("f0")
vds = hl.vds.filter_samples(vds, samples_table, keep=True)

# Filter variants to interval
vds = hl.vds.filter_variants(vds, lambda locus: interval.contains(locus))

# Densify only now
mt = hl.vds.to_dense_mt(vds)

# Filter to PASS variants
mt = mt.filter_rows(mt.filters.length() == 0)

# Export to VCF
hl.export_vcf(mt, output_path)
CODE
  >>>

  runtime {
    docker: "hailgenetics/hail:0.2.134-py3.11"
    memory: "16G"
    cpu: 2
    disks: "local-disk 100 HDD"
  }

  output {
    File vcf_output = "${output_bucket}/region_${interval.replace(':', '_').replace('-', '_')}.vcf.bgz"
  }
}
