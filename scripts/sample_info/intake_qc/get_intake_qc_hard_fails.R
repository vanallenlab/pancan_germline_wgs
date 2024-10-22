#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Compile list of samples IDs that hard fail intake QC
# The reasons for hard failure are:
# 1. inferred sex != reported sex
# 2. ambiguous inferred sex
# 3. complete lack of any phenotype data


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
parser <- ArgumentParser(description="Find intake QC hard fail samples")
parser$add_argument("--qc-tsv", metavar=".tsv", type="character",
                    help="Intake QC .tsv", required=TRUE)
parser$add_argument("--outfile", metavar=".txt", type="character", required=TRUE,
                    help="Path to output .txt")
args <- parser$parse_args()

# # DEV:
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.post_qc_batching.tsv.gz",
#              "outfile" = "~/scratch/dfci-g2c.intake_qc.hard_fail.samples.list")

# Load QC matrix
qc.df <- read.table(args$qc_tsv, header=T, sep="\t", comment.char="", check.names=F)

# Find samples with inferred & reported sex mismatch
sex.mm.fails <- qc.df[which(qc.df$batching_sex != qc.df$reported_sex
                            & qc.df$sex_karyotype %in% c("XX", "XY")), 1]
cat(paste(prettyNum(length(sex.mm.fails), big.mark=","), "samples had mismatching",
          "inferred and reported sexes\n"))

# Find samples with ambiguous sexes
sex.amb.margin <- 0.15
sex.amb.fails <- qc.df[which(qc.df$chrX_ploidy >= 1 + sex.amb.margin
                             & qc.df$chrX_ploidy <= 2 - sex.amb.margin
                             & qc.df$chrY_ploidy >= sex.amb.margin
                             & qc.df$chrY_ploidy <= 1 - sex.amb.margin), 1]
cat(paste(prettyNum(length(sex.amb.fails), big.mark=","), "samples had ambiguous",
          "inferred sexes\n"))

# Find samples completely lacking all phenotype data
pheno.cols <- c("reported_sex", "cancer", "age", "years_to_last_contact",
                "years_left_censored", "vital_status", "stage")
no.pheno.fails <- qc.df[which(apply(qc.df[, pheno.cols], 1, function(v){all(is.na(v))})), 1]
cat(paste(prettyNum(length(no.pheno.fails), big.mark=","), "samples had a complete",
          "lack of any phenotype data available\n"))

# Unify three causes of hard failures and write to --outfile
all.fails <- sort(unique(c(sex.mm.fails, sex.amb.fails, no.pheno.fails)))
cat(paste(prettyNum(length(all.fails), big.mark=","), "samples failed for",
          "any of the above three reasons\n"))
write.table(all.fails, args$outfile, col.names=F, row.names=F, quote=F)
