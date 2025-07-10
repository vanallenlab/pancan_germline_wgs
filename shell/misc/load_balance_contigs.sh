#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Load-balance sets of chromosomes for optimal parallelization across workspaces

# Resolve the absolute path to this script
# Only necessary when executed from the command line
# When run interactively, you need to set SCRIPT_DIR to the path of the directory containing this script
SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do
  DIR="$(cd -P "$(dirname "$SOURCE")" >/dev/null 2>&1 && pwd)"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE"
done
SCRIPT_DIR="$(cd -P "$(dirname "$SOURCE")" >/dev/null 2>&1 && pwd)"

# Build local execution directory
WRKDIR=`mktemp -d`
cd $WRKDIR >/dev/null 2>&1

# Localize files referenced multiple times
gsutil cp gs://dfci-g2c-refs/hg38/hg38.genome ./

# Compute four statistics for each contig:
# 1. Total alignable length in hg38 primary assembly
wget -nv -O - \
  https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz \
| gunzip -c | cut -f2-4 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - \
| bgzip -c > hg38.gap.bed.gz
awk -v OFS="\t" '{ print $1, "0", $2 }' hg38.genome | sort -Vk1,1 \
| bedtools subtract -a - -b hg38.gap.bed.gz \
| awk -v FS="\t" -v OFS="\t" \
  '{ sum[$1] += $3 - $2 } END { for (chr in sum) print chr, sum[chr] }' \
| sort -k1,1 \
> hg38.alignable_bp.tsv

# 2. Total number of protein-coding genes
wget -nv -O - \
  https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.4/MANE.GRCh38.v1.4.ensembl_genomic.gtf.gz \
| gunzip -c | fgrep -w "gene_type \"protein_coding\"" \
| fgrep -wf <( cut -f1 hg38.genome ) \
| awk -v FS="\t" -v OFS="\t" \
  '{ if ($3=="gene") count[$1] += 1 } END { for (chr in count) print chr, count[chr] }' \
| sort -k1,1 \
> hg38.gene_count.tsv

# 3. Total number of SNVs/indels in gnomAD
while read contig; do
  gsutil -m cat \
    gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/$contig/gnomad.v4.1.$contig.af_distribution.merged.tsv.gz \
  | gunzip -c | grep -ve '^#' \
  | awk -v FS="\t" -v OFS="\t" -v contig=$contig \
    '$1 != "sv" { for (i = 3; i <= NF; i++) sum += $i } END { print contig, sum }'
done < <( cut -f1 hg38.genome | sort -k1,1 ) \
> hg38.snv_count.tsv

# 4. Total number of SVs in gnomAD
while read contig; do
  gsutil -m cat \
    gs://dfci-g2c-refs/gnomad/gnomad_v4_site_metrics/$contig/gnomad.v4.1.$contig.af_distribution.merged.tsv.gz \
  | gunzip -c | grep -ve '^#' \
  | awk -v FS="\t" -v OFS="\t" -v contig=$contig \
    '$1 == "sv" { for (i = 3; i <= NF; i++) sum += $i } END { print contig, sum }'
done < <( cut -f1 hg38.genome | sort -k1,1 ) \
> hg38.sv_count.tsv

# Merge metrics collected above
join -j 1 -t$'\t' hg38.alignable_bp.tsv hg38.gene_count.tsv \
| join -j 1 -t$'\t' - hg38.snv_count.tsv \
| join -j 1 -t$'\t' - hg38.sv_count.tsv \
| sort -Vk1,1 \
| cat <( echo -e "#contig\tbp\tgenes\tsnvs\tsvs" ) - \
| gzip -c \
> hg38.load_balancing_metrics.tsv.gz

# Divide contigs into five groups that are approximately balanced across all four metrics above
$SCRIPT_DIR/../../scripts/utilities/loadBalancedSplitter.R \
  --input-tsv hg38.load_balancing_metrics.tsv.gz \
  --number-of-splits 5 \
  --out-prefix hg38.split
for i in $( seq 1 5 ); do
  mv hg38.split$i dfci-g2c.v1.contigs.w$i.list
done

# Write chr19, chr21, chr22, and chrY to a small "development" list
# This is only for development purposes of GATK-HC joint genotyping and downstream QC
cat << EOF > dfci-g2c.v1.contigs.dev.list
chr19
chr21
chr22
chrY
EOF

# Copy all contig lists to central directory
gsutil -m cp \
  dfci-g2c.v1.contigs.*.list \
  gs://dfci-g2c-refs/hg38/contig_lists/

# Copy load balancing metrics to central directory
gsutil -m cp \
  hg38.load_balancing_metrics.tsv.gz \
  gs://dfci-g2c-refs/hg38/

# Switch back to execution directory and delete working directory
# (This is only used if/when this script is executed from the command line)
cd - >/dev/null 2>&1
rm -rf $WRKDIR
