#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Update G2C sample metadata with a new QC failure column


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(RLCtools, quietly=TRUE)


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Append a QC failure column to G2C sample metadata")
parser$add_argument("--qc-tsv", metavar=".tsv", type="character",
                    help="Intake QC .tsv", required=TRUE)
parser$add_argument("--new-column-name", metavar="string", type="character",
                    help="Name of new column to be added to metadata", required=TRUE)
parser$add_argument("--pass-samples-list", metavar=".txt", type="character",
                    help="List of samples passing QC")
parser$add_argument("--fail-samples-list", metavar=".txt", type="character",
                    help="List of samples failing QC")
parser$add_argument("--all-samples-list", metavar=".txt", type="character",
                    help="List of samples considered during QC")
parser$add_argument("--outfile", metavar=".tsv", type="character", required=TRUE,
                    help="Path to output .tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.post_qc_batching.tsv.gz",
#              "pass_samples_list" = "~/scratch/dev.qc.pass.samples.list",
#              "new_column_name" = "new_qc_pass",
#              "outfile" = "~/scratch/dev.qc_update.w_new_fail.tsv")

# One of --pass-samples-list or --fail-samples-list must be provided
if(is.null(args$pass_samples_list) & is.null(args$fail_samples_list)){
  stop("Must provide either --pass-samples-list or --fail-samples-list")
}

# Load QC matrix
qc.df <- read.table(args$qc_tsv, header=T, sep="\t", comment.char="", check.names=F)

# Determine universe of all eligible samples
all.ids <- if(!is.null(args$all_samples_list)){
  sort(unique(read.table(args$all_samples_list, header=F)[, 1]))
}else{sort(unique(qc.df[, 1]))}

# Load list of samples passing QC
pass.ids <- if(!is.null(args$pass_samples_list)){
  sort(unique(read.table(args$pass_samples_list, header=F)[, 1]))
}else{c()}

# Load list of samples failing QC
fail.ids <- if(!is.null(args$fail_samples_list)){
  sort(unique(read.table(args$fail_samples_list, header=F)[, 1]))
}else{c()}

# Check that all pass/fail samples are present in universe
if(length(setdiff(c(pass.ids, fail.ids), all.ids)) > 0){
  stop(paste("Not all samples listed in either --pass-samples-list and/or",
             "--fail-samples-list are present in --all-samples-list",
             "and/or --qc-tsv"))
}

# Ensure that no samples are specified as both passing and failing
if(length(intersect(pass.ids, fail.ids)) > 0){
  stop(paste("Sample overlaps between --pass-samples-list and",
             "--fail-samples-list are not permitted"))
}

# Annotate metadata with new QC pass column
# Using pythonic boolean strings for consistency with prior QC
qc.df[, args$new_column_name] <- NA
qc.df[which(qc.df[, 1] %in% all.ids), args$new_column_name] <- "True"
qc.df[which(qc.df[, 1] %in% fail.ids), args$new_column_name] <- "False"
qc.df[which(qc.df[, 1] %in% pass.ids), args$new_column_name] <- "True"

# Write new metadata to --outfile
write.table(qc.df, args$outfile, col.names=T, row.names=F, sep="\t", quote=F)
