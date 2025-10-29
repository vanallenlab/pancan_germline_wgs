#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Divide rows from a .tsv matrix into groups by load balancing features per row


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)


#############
# Functions #
#############
# Divide a normalized matrix into feature-balanced groups
load.balance <- function(x, n.splits){

  # Record original order for sorting outputs
  orig.order <- rownames(x)

  # Reorder matrix by largest minimum feature
  x <- x[order(apply(x, 1, min), decreasing=TRUE), ]

  # Initialize split collector with one of the largest elements per split
  splits <- lapply(1:n.splits, function(i){rownames(x)[i]})
  x.remain <- x[-c(1:n.splits), ]

  # Iterate over remaining rows
  while(nrow(x.remain) > 0){
    next.key <- rownames(x.remain)[1]
    s.totals <- lapply(splits, function(k){apply(x[k, ], 2, sum)})
    add.to <- which.min(sapply(s.totals, function(v){
      min(as.numeric(v) + as.numeric(x[next.key, ]))
    }))
    splits[[add.to]] <- c(splits[[add.to]], next.key)
    x.remain <- x.remain[-c(1), ]
  }

  # Return optimized splits, sorted according to input order
  lapply(splits, function(v){
    names(sort(sapply(v, function(vi){which(orig.order == vi)})))
  })
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Partition elements with feature load balancing")
parser$add_argument("--input-tsv", required=TRUE, metavar=".tsv",
                    help=paste("Input .tsv to be partitioned. First column must",
                               "be unique key. All other columns must be numeric."))
parser$add_argument("-n", "--number-of-splits", metavar="integer",
                    type="numeric", required=TRUE,
                    help=paste("Number of groups to create"))
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output .txt files",
                    default="./load_balanced_split")
args <- parser$parse_args()

# DEV:
# args <- list("input_tsv" = "~/scratch/hg38.load_balancing_metrics.tsv.gz",
#              "number_of_splits" = 5,
#              "out_prefix" = "~/scratch/hg38.load_balanced.group")

# Read data
x <- read.table(args$input_tsv, sep="\t", header=F)
n.features <- ncol(x)-1
if(n.features == 0){stop("No features found in --input-tsv")}
x[, -c(1)] <- apply(x[, -c(1)], 2, as.numeric)
rownames(x) <- x[, 1]
x[, 1] <- NULL
colnames(x) <- paste("f", 1:n.features, sep="")
n.obs <- nrow(x)

# Proceed based on number of desired splits
if(args$number_of_splits == 0){

  stop("Must specify a positive integer for --number-of-splits")

}else if(args$number_of_splits == 1){

  splits <- list(rownames(x))

}else if(args$number_of_splits >= n.obs){

  warning(paste("--number-of-splits exceeds number of rows in --input-tsv;",
                "writing each element to its own file"))
  splits <- as.list(rownames(x))

}else{

  # Normalize all features
  x[, 1:n.features] <- apply(x[, 1:n.features], 2, function(v){v/sum(v)})

  # Load balance groups across all features
  splits <- load.balance(x, args$number_of_splits)

}

# Write splits to output file(s)
sink <- sapply(1:length(splits), function(s){
  write.table(splits[[s]], paste(args$out_prefix, s, sep=""),
              col.names=F, row.names=F, quote=F)
})
