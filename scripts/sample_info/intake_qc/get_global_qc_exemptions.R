#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Compile list of samples IDs that are exempt from certain global QC metrics
# These include:
# 1. HMF samples are exempt from WGD
# 2. Samples with exactly one unambiguous aneuploidy are exempt from
#    SNV allele balance-based metrics


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(RLCtools, quietly=TRUE)


##################
# Data Functions #
##################
# Get list of sample IDs with extremely clear & isolated autosomal aneuploidies
get.autosomal.aneuploidies <- function(qc.df, tolerance=0.25){
  auto.ploidy.cols <- paste("chr", 1:22, "_ploidy", sep="")
  has.clean.aneu <- apply(qc.df[, auto.ploidy.cols], 1, function(v){
    ploidy.resid <- abs(v - 2)
    (max(ploidy.resid) >= 1-tolerance
     & max(ploidy.resid[-which.max(ploidy.resid)]) <= tolerance)
  })
  qc.df$G2C_id[which(has.clean.aneu)]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Find global QC exemptions")
parser$add_argument("--qc-tsv", metavar=".tsv", type="character",
                    help="Intake QC .tsv", required=TRUE)
parser$add_argument("--outfile", metavar=".tsv", type="character", required=TRUE,
                    help="Path to output .tsv of exemptions")
args <- parser$parse_args()

# # DEV:
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.post_qc_batching.tsv.gz",
#              "outfile" = "~/scratch/dfci-g2c.intake_qc.global_exemptions.tsv")

# Load QC matrix
qc.df <- read.table(args$qc_tsv, header=T, sep="\t", comment.char="", check.names=F)
colnames(qc.df) <- gsub("#", "", colnames(qc.df), fixed=T)

# 1. All HMF samples are exempt from WGD
hmf.ids <- qc.df$G2C_id[which(qc.df$cohort == "hmf")]
ex.df <- data.frame("sample"=hmf.ids, "metric"=rep("wgd_score", length(hmf.ids)))
cat(paste(prettyNum(length(hmf.ids), big.mark=","), "samples from HMF are being",
          "exempted from global WGD-based QC\n"))

# 2a. Samples with sex chromosome aneuploidies and no autosomal aneuploidies
#     are exempt from allele balance-based QC
ab.qc <- c("charr", "mean_ref_ab_hom_alt", "inconsistent_ab_het_rate")
sex.aneu.ids <- qc.df$G2C_id[which(!qc.df$sex_karyotype %in% c("XX", "XY"))]
sex.aneu.df <- data.frame("sample"=stretch.vector(sex.aneu.ids, length(ab.qc)),
                          "metric"=rep(ab.qc, length(sex.aneu.ids)))
ex.df <- rbind(ex.df, sex.aneu.df)
cat(paste(prettyNum(length(sex.aneu.ids), big.mark=","), "samples with sex",
          "chromosome aneuploidies are being exempted from SNV allele",
          "balance-based QC\n"))

# 2b. Samples with very high-quality autosomal aneuploidies are exempt from
#     are exempt from allele balance-based QC
auto.aneu.ids <- get.autosomal.aneuploidies(qc.df)
auto.aneu.df <- data.frame("sample"=stretch.vector(auto.aneu.ids, length(ab.qc)),
                          "metric"=rep(ab.qc, length(auto.aneu.ids)))
ex.df <- rbind(ex.df, auto.aneu.df)
cat(paste(prettyNum(length(auto.aneu.ids), big.mark=","), "samples with clean",
          "autosomal aneuploidies are being exempted from SNV allele",
          "balance-based QC\n"))

# Unify three causes of hard failures and write to --outfile
write.table(ex.df, args$outfile, sep="\t", col.names=F, row.names=F, quote=F)
