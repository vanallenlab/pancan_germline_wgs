#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Filter and sample a fixed number of samples from a priority list for QC


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
parser <- ArgumentParser(description="Select samples to use for sample-level QC")
parser$add_argument("--all-samples-list", metavar=".txt", type="character", required=TRUE,
                    help="List of all eligible sample IDs")
parser$add_argument("-N", metavar="int", type="numeric", required=TRUE,
                    help="Number of samples to select")
parser$add_argument("--priority-tsv", metavar=".tsv", type="character",
                    help=paste("Optional three-column .tsv of sample IDs, ",
                               "priority tier, and sampling weight. Samples ",
                               "present in --all-samples-list but missing from",
                               "--priority-tsv will be imputed with tier=0 and",
                               "weight=1"))
parser$add_argument("--out-list", metavar=".txt", type="character", required=TRUE,
                    help="path to output list of selected samples")
args <- parser$parse_args()

# # DEV (AoU RW)
# args <- list("all_samples_list" = "staging/sample_priority/g2c_ids.present_after_calling.samples.list",
#              "N" = 2000,
#              "priority_tsv" = "staging/sample_priority/dfci-g2c.v1.sample_qc_priority.tsv",
#              "out_list" = "scratch/n2k.test.samples.list")

# Read list of all samples
s.all <- sort(unique(read.table(args$all_samples_list, header=F)[, 1]))

# Import sample priority, if optioned
# Otherwise, set all samples to equal priority
if(!is.null(args$priority_tsv)){
  s.p <- read.table(args$priority_tsv, header=F, sep="\t")
  colnames(s.p) <- c("sid", "tier", "weight")
  s.p <- s.p[which(s.p$sid %in% s.all), ]
  sid.missing <- s.all[which(!(s.all %in% s.p$sid))]
  if(length(sid.missing) > 0){
    s.p.add <- data.frame("sid"=sid.missing, "tier"=0, "weight"=1)
    s.p <- as.data.frame(rbind(s.p, s.p.add))
  }
}else{
  s.p <- data.frame("sid"=s.all, "tier"=0, "weight"=1)
}
s.p <- s.p[with(s.p, order(-tier, -weight, sid)), ]

# Sample desired number of IDs
n.per.tier <- table(s.p$tier)
n.per.tier <- n.per.tier[order(-as.numeric(names(n.per.tier)))]
auto.tiers <- names(n.per.tier)[which(cumsum(n.per.tier) <= args$N)]
if(length(auto.tiers) > 0){
  s.keep <- s.p[which(s.p$tier %in% auto.tiers), "sid"]
}else{
  s.keep <- c()
}
random.tier <- names(n.per.tier)[which.min(cumsum(n.per.tier) <= args$N)]
n.to.sample <- args$N - length(s.keep)
set.seed(2025)
s.add <- sample(s.p[which(s.p$tier == random.tier), "sid"], n.to.sample,
                replace=FALSE, prob=s.p[which(s.p$tier == random.tier), "weight"])
s.keep <- sort(unique(c(s.keep, s.add)))

# Write list of sampled IDs to --out-list
write.table(s.keep, args$out_list, col.names=F, row.names=F, quote=F)
