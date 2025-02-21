#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Update phenotype data in the overall G2C sample metadata .tsv


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Update sample phenotypes in G2C metadata.tsv")
parser$add_argument("--metadata-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Main G2C sample metadata .tsv for all samples")
parser$add_argument("--phenotypes-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Sample phenotypes .tsv to use for updates")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("metadata_tsv" = "~/scratch/dfci-g2c.intake_qc.all.post_qc_batching.non_aou.tsv.gz",
#              "phenotypes_tsv" = "~/scratch/ceph.phenos.tsv.gz",
#              "out_tsv" = "~/scratch/dfci-g2c.intake_qc.all.post_qc_batching.non_aou.ceph_updated.tsv")

# Load sample metadata, leaving all column names untouched
main.df <- read.table(args$metadata_tsv, header=T, sep="\t",
                      comment.char="", check.names=F)

# Load phenotypes and adjust columns for merging
update.df <- read.table(args$phenotypes_tsv, sep="\t", quote="",
                        comment.char="", check.names=F, header=T)
colnames(update.df)[which(colnames(update.df) == "Sample")] <- "original_id"
colnames(update.df)[which(colnames(update.df) == "Cohort")] <- "cohort"
shared.columns <- intersect(colnames(main.df), colnames(update.df))
update.df <- update.df[, shared.columns]

# Add helper index to both dataframes for merging
rownames(main.df) <- apply(main.df[, c("cohort", "original_id")],
                           1, paste, collapse="|")
rownames(update.df) <- apply(update.df[, c("cohort", "original_id")],
                             1, paste, collapse="|")

# Update main dataframe
update.ids <- intersect(rownames(main.df), rownames(update.df))
main.df[update.ids, shared.columns] <- update.df[update.ids, shared.columns]

# Write updated manifest to --out-tsv
write.table(main.df, args$out_tsv, col.names=T, row.names=F, quote=F, sep="\t")
