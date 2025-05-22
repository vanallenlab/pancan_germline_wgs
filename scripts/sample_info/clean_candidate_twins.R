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
twins <- read.table(args[1], header=T, sep="\t", comment.char="")
