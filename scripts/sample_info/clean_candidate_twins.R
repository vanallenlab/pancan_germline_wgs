#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Remove likely duplicated files from a list of candidate twins/replicates

options(scipen=1000, stringsAsFactors=F)

# Three positional arguments:
# 1. Input .tsv of candidate twins from InferTwins.wdl
# 2. Input .tsv of sample metadata
# 3. Path to output .tsv of filtered twins file
args <- commandArgs(trailingOnly=TRUE)

# Read inputs
twins <- read.table(args[1], header=T, sep="\t", comment.char="", check.names=F)
meta <- read.table(args[2], header=T, sep="\t", comment.char="", check.names=F)
rownames(meta) <- meta$G2C_id
meta$G2C_id <- NULL

# Filter to twins with kinship >=0.45
twins <- twins[which(twins$KINSHIP >= 0.45), ]

# Subset metadata to technical values and standard normalize
qc.cols <- c("grafpop_GD1", "grafpop_GD2", "grafpop_GD3", "chrX_ploidy",
             "chrY_ploidy", "charr", "insert_size", "read_length",
             "mean_coverage", "median_coverage", "wgd_score")
meta[, qc.cols] <- apply(meta[, qc.cols], 2, scale, center=TRUE, scale=TRUE)
meta <- meta[, qc.cols]

# Compute technical similarity for all pairs of twins
qc.d <- sapply(1:nrow(twins), function(i){
  tids <- as.character(twins[i, c("IID1", "IID2")])
  v1 <- as.numeric(meta[tids[1], ])
  v2 <- as.numeric(meta[tids[2], ])
  sqrt(sum((v1 - v2)^2))
})

# Only keep twins with Euclidean distance of QC Z-scores > 0.5
write.table(twins[which(qc.d > 0.5), ], args[3], col.names=T, row.names=F,
            quote=F, sep="\t")
