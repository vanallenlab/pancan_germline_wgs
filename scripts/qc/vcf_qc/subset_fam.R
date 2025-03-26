#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Subset a .fam file to duos or trios


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
parser <- ArgumentParser(description="Subset a .fam file to duos or trios")
parser$add_argument("--in-fam", metavar=".fam", type="character", required=TRUE,
                    help="Input .fam file")
parser$add_argument("--all-samples", metavar="string", type="character",
                    help="Optional list of samples to include when filtering")
parser$add_argument("--input-has-header", action="store_true", default=FALSE,
                    help="Input .fam file has a header line")
parser$add_argument("--out-fam", metavar=".fam", type="character", required=TRUE,
                    help="path to output .fam")
args <- parser$parse_args()

# # DEV:
# args <- list("in_fam" = "~/scratch/dfci-g2c.hgsvc.fam",
#              "all_samples" = "~/scratch/test.pass.samples.list",
#              "input_has_header" = FALSE,
#              "out_fam" = "~/scratch/dfci-g2c.hgsvc.subset_test.fam")

# Load input .fam
fam <- G2CR::load.famfile(args$in_fam, header=args$input_has_header)

# Apply --all-samples restriction, if optioned
if(!is.null(args$all_samples)){
  pass.ids <- unique(read.table(args$all_samples, header=F)[, 1])
  names(pass.ids) <- pass.ids
  id.cols <- c("proband", "father", "mother")
  fam[, id.cols] <- apply(fam[, id.cols], 2, RLCtools::remap, map=pass.ids, default.value=0)
}

# Restrict to duos and trios
has.pro <- !(fam$proband %in% c(0, ".") | is.na(fam$proband))
has.fa <- !(fam$father %in% c(0, ".") | is.na(fam$father))
has.mo <- !(fam$mother %in% c(0, ".") | is.na(fam$mother))
fam <- fam[which(has.pro & (has.fa | has.mo)), ]

# Write updated manifest to --out-fam
write.table(fam, args$out_fam, col.names=F, row.names=F, quote=F, sep="\t")
