#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Merge sub-components of single-cohort intake QC metrics
# Helper script for consolidate_cohort_intake_qc.sh


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Merge intake QC sub-components for a single cohort")
parser$add_argument("--demo", metavar=".tsv", type="character", required=TRUE,
                    help="Sample demographics")
parser$add_argument("--charr", metavar=".tsv", type="character", required=TRUE,
                    help="Charr metrics")
parser$add_argument("--cov", metavar=".tsv", type="character", required=TRUE,
                    help="GATK-SV coverage means & medians")
parser$add_argument("--reads", metavar=".tsv", type="character", required=TRUE,
                    help="GATK read lengths & insert sizes")
parser$add_argument("--ploidy", metavar=".tsv", type="character", required=TRUE,
                    help="GATK-SV ploidy & WGD")
parser$add_argument("--outfile", metavar="path", type="character", required=TRUE,
                    help="Path to output .tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("demo" = "/var/folders/zj/59lpkx5926376cnm6l48d2x00000gp/T/apollo_demo/apollo.demo.tsv.gz",
#              "charr" = "/var/folders/zj/59lpkx5926376cnm6l48d2x00000gp/T/apollo_charr/apollo.charr.tsv.gz",
#              "cov" = "/var/folders/zj/59lpkx5926376cnm6l48d2x00000gp/T/apollo_cov/apollo.cov.tsv.gz",
#              "reads" = "/var/folders/zj/59lpkx5926376cnm6l48d2x00000gp/T/apollo_read/apollo.read_metrics.tsv.gz",
#              "ploidy" = "/var/folders/zj/59lpkx5926376cnm6l48d2x00000gp/T/apollo.ploidy.tsv.gz",
#              "outfile" = "~/scratch/apollo.intake_qc.test.tsv")

# Load & clean demographic data
demo.df <- read.table(args$demo, check.names=F, header=T, sep="\t", comment.char="")
demo.col.rename <- c("SNPs" = "grafpop_SNPs",
                     "GD1(x)" = "grafpop_GD1",
                     "GD2(y)" = "grafpop_GD2",
                     "GD3(z)" = "grafpop_GD3",
                     "E(%)" = "pct_EUR",
                     "F(%)" = "pct_AFR",
                     "A(%)" = "pct_ASN")
for(i in 1:length(demo.col.rename)){
  old.name <- names(demo.col.rename)[i]
  new.name <- demo.col.rename[i]
  colnames(demo.df)[which(colnames(demo.df) == old.name)] <- new.name
}

# Load other data
charr.df <- read.table(args$charr, check.names=F, header=T, sep="\t", comment.char="")
cov.df <- read.table(args$cov, check.names=F, header=T, sep="\t", comment.char="")
reads.df <- read.table(args$reads, check.names=F, header=T, sep="\t", comment.char="")
reads.df$insert_size <- round(as.numeric(reads.df$insert_size))
colnames(reads.df)[1] <- "Sample"
ploidy.df <- read.table(args$ploidy, check.names=F, header=T, sep="\t", comment.char="")
colnames(ploidy.df)[1:2] <- c("Cohort", "Sample")

# Merge data
out.df <- merge(demo.df, charr.df, by.x="Sample", by.y="#SAMPLE", all=F, sort=F)
out.df <- merge(out.df, cov.df, by="Sample", all.x=T, all.y=F, sort=F)
out.df <- merge(out.df, reads.df, by="Sample", all.x=T, all.y=F, sort=F)
out.df <- merge(out.df, ploidy.df, by="Sample", all=F, sort=F)

# Write to --out-tsv
out.col.order <- c("Sample", "Cohort", "grafpop_SNPs", "grafpop_GD1", "grafpop_GD2",
                   "grafpop_GD3", "pct_EUR", "pct_AFR", "pct_ASN", "chrX_count",
                   "chrY_count", "ancestry", "HQ_HOM", "HQ_HOM_RATE", "HQ_HET",
                   "HQ_HET_RATE", "CHARR", "MEAN_REF_AB_HOM_ALT", "HETEROZYGOSITY_RATE",
                   "INCONSISTENT_AB_HET_RATE", "insert_size", "read_length",
                   "rd_median", "rd_mean", "chr1_CopyNumber", "chr2_CopyNumber",
                   "chr3_CopyNumber", "chr4_CopyNumber", "chr5_CopyNumber",
                   "chr6_CopyNumber", "chr7_CopyNumber", "chr8_CopyNumber",
                   "chr9_CopyNumber", "chr10_CopyNumber", "chr11_CopyNumber",
                   "chr12_CopyNumber", "chr13_CopyNumber", "chr14_CopyNumber",
                   "chr15_CopyNumber", "chr16_CopyNumber", "chr17_CopyNumber",
                   "chr18_CopyNumber", "chr19_CopyNumber", "chr20_CopyNumber",
                   "chr21_CopyNumber", "chr22_CopyNumber", "chrX_CopyNumber",
                   "chrY_CopyNumber", "median_coverage", "wgd_score", "nondiploid_bins")
write.table(out.df[, out.col.order], args$outfile, col.names=T, row.names=F,
            sep="\t", quote=F)
