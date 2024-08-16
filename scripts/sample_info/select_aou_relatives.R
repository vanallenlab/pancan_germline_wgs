#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Select relatives of G2C cancer cases from AoU to be included in phase 1


# Set parameters & load libraries
options(scipen=1000, stringsAsFactors=FALSE)
require(argparse, quietly=TRUE)

# Parse command line arguments and options
parser <- ArgumentParser(description="Select relatives from All of Us")
parser$add_argument("--kinship", metavar=".tsv", type="character",
                    help="Precompued AoU relatedness.tsv", required=TRUE)
parser$add_argument("--case-ids", metavar=".txt", type="character",
                    help="List of AoU sample IDs for G2C cases", required=TRUE)
parser$add_argument("--ufc-case-ids", metavar=".txt", type="character",
                    help="List of AoU sample IDs for UFC cases", required=TRUE)
parser$add_argument("--processed-ids", metavar=".txt", type="character",
                    help="List of samples that have been processed", required=TRUE)
parser$add_argument("--max-new-samples", metavar="integer", default=400,
                    help="Maximum number of new samples to process.")
parser$add_argument("--outfile", metavar="path", type="character",
                    help="path to output list of relative sample IDs", required=TRUE)
args <- parser$parse_args()

# # DEV (only works on AoU RW):
# args <- list("kinship" = "relatedness.tsv",
#              "case_ids" = "all_g2c_cases.samples.list",
#              "ufc_case_ids" = "ufc_cohort.list",
#              "processed_ids" = "all_g2c.samples.list",
#              "max_new_samples" = 400,
#              "outfile" = "sample_lists/relatives.samples.list")

# Read kinship data
pairs <- read.table(args$kinship, header=T, sep="\t")
colnames(pairs) <- c("s1", "s2", "kin")

# Annotate kinship data with G2C case overlap
case.ids <- unique(read.table(args$case_ids, header=F)[, 1])
pairs$n_case <- apply(pairs[, 1:2], 1, function(sids){length(intersect(sids, case.ids))})

# Annotate kinship data with UFC case overlap
ufc.case.ids <- unique(read.table(args$ufc_case_ids, header=F)[, 1])
pairs$n_ufc_case <- apply(pairs[, 1:2], 1, function(sids){length(intersect(sids, ufc.case.ids))})

# Annotate kinship data with number of samples already processed for G2C
processed.ids <- unique(read.table(args$processed_ids, header=F)[, 1])
pairs$n_processed <- apply(pairs[, 1:2], 1, function(sids){length(intersect(sids, processed.ids))})

# Count number of first- and second-degree relatives per case
cases.in.pairs <- case.ids[which(case.ids %in% unlist(pairs[, 1:2]))]
fdr.per.case <- sapply(cases.in.pairs, function(sid){
  max(c(0, length(unique(unlist(pairs[which(pairs$kin > 0.2 & (pairs$s1 == sid | pairs$s2 == sid)), 1:2]))) - 1))
})
names(fdr.per.case) <- cases.in.pairs
sdr.per.case <- sapply(cases.in.pairs, function(sid){
  max(c(0, length(unique(unlist(pairs[which(pairs$kin <= 0.2 & (pairs$s1 == sid | pairs$s2 == sid)), 1:2]))) - 1))
})
names(sdr.per.case) <- cases.in.pairs

# Annotate kinship data with total number of relatives of cases involved in each pair
pairs$n_fdr <- apply(pairs[, 1:2], 1, function(sids){
  sum(c(0, fdr.per.case[as.character(intersect(sids, cases.in.pairs))]))
})
pairs$n_sdr <- apply(pairs[, 1:2], 1, function(sids){
  sum(c(0, sdr.per.case[as.character(intersect(sids, cases.in.pairs))]))
})

# Randomly shuffle pairs before sorting by prioritization criteria
set.seed(2024)
pairs <- pairs[sample(1:nrow(pairs), nrow(pairs), replace=F), ]
pairs <- pairs[with(pairs, order(-n_case, -n_ufc_case, -n_fdr, -n_sdr, -n_processed, -kin)), ]

# Iterate over priority-ordered pairs, adding all relatives of any cases in each row
remaining.pairs <- pairs
all_noncase <- all_new <- c()
n_noncase <- n_new <- 0
while(n_new < args$max_new_samples){
  next.cases <- intersect(as.character(unlist(remaining.pairs[1, 1:2])), case.ids)
  hit.rels <- setdiff(unlist(remaining.pairs[which(remaining.pairs$s1 %in% next.cases
                                                   | remaining.pairs$s2 %in% next.cases),
                                             1:2]), case.ids)
  all_noncase <- unique(c(all_noncase, hit.rels))
  n_noncase <- length(all_noncase)
  all_new <- setdiff(all_noncase, processed.ids)
  n_new <- length(all_new)
  remaining.pairs <- remaining.pairs[-which(remaining.pairs$s1 %in% next.cases
                                           | remaining.pairs$s2 %in% next.cases), ]
}

# Count number of unique cases with at least one non-case relative to be included
cases_w_rels <- intersect(unlist(pairs[which(pairs$s1 %in% all_noncase
                                             | pairs$s2 %in% all_noncase),
                                       1:2]),
          case.ids)
n_cases_w_rels <- length(cases_w_rels)

# Print summary results
cat("Relative selection complete. Results:\n")
cat(paste("Total number of non-case relatives included:", n_noncase, "\n"))
cat(paste("Number of non-case relatives to be processed:", n_new, "\n"))
cat(paste("Number of cases with at least one relative included:", n_cases_w_rels, "\n"))

# Write list of relatives to be processed to --outfile
write.table(sort(all_new), args$outfile, quote=F, row.names=F, col.names=F)

