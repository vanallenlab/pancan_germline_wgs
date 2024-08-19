#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Merge intake QC metrics with basic phenotype data required for QC & batching


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2C, quietly=TRUE)


##################
# Data Functions #
##################
# Load phenotype data and subset only to columns needed for QC and batching
load.pheno.df <- function(tsv.in){
  df <- read.table(tsv.in, header=T, sep="\t")
  keep.cols <- c("Sample", "Cohort", "reported_sex", "cancer")
  df[, keep.cols]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Visualize intake QC metrics")
parser$add_argument("--qc-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Intake QC .tsv")
parser$add_argument("--phenos-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Sample phenotype .tsv")
parser$add_argument("--out-tsv", metavar="path", type="character", required=TRUE,
                    help="Path to output .tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.tsv.gz",
#              "phenos_tsv" = "~/scratch/dfci-g2c.non_aou.phenos.tsv.gz",
#              "out_tsv" = "~/scratch/dfci-g2c.intake_qc.local_test.wphenos.tsv")

# Load data
qc.df <- read.table(args$qc_tsv, check.names=F, header=T, sep="\t", comment.char="")
pheno.df <- load.pheno.df(args$phenos_tsv)

# Merge data
out.df <- merge(qc.df, pheno.df, by=c("Sample", "Cohort"), all.x=T, all.y=F, sort=F)

# Write to --out-tsv
write.table(out.df, args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
