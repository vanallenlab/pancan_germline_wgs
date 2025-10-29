#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Produce a .tsv mapping G2C IDs to external IDs for a single cohort


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)


##################
# Data functions #
##################
# Load a .kin0 file of likely twins
load.kin0.twins <- function(twins.in){
  read.table(twins.in, header=T, sep="\t", comment.char="")[, c(2, 4)]
}

# "Explode" ID map for pairs of inter-cohort twins
explode.twins <- function(df, twin.map){
  add.rows <- do.call("rbind", lapply(1:nrow(twin.map), function(i){
    tids <- as.character(unlist(twin.map[i, ]))
    hit.idx <- tids %in% df[, 1]
    if(sum(hit.idx) == 1){
      anchor.id <- tids[which(hit.idx)]
      satellite.id <- tids[which(!hit.idx)]
      external.id <- df[which(df[, 1] == anchor.id), 2]
      c(satellite.id, external.id)
    }
  }))
  if(nrow(add.rows) > 0){
    add.rows <- as.data.frame(add.rows)
    colnames(add.rows) <- colnames(df)
    as.data.frame(rbind(df, add.rows))
  }else{
    return(df)
  }
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Map G2C IDs to external IDs")
parser$add_argument("--metadata-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Main G2C sample metadata .tsv for all samples")
parser$add_argument("--cohort", type="character", metavar="string", required=TRUE,
                    help="Cohort of interest")
parser$add_argument("--twins", metavar=".kin0", help=".kin0 file of putative twins")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to id_map.tsv")
args <- parser$parse_args()

# # DEV (AoU RW):
# args <- list("metadata_tsv" = "staging/external_id_maps/dfci-g2c.sample_meta.gatkhc_posthoc_outliers.tsv.gz",
#              "cohort" = "hgsvc",
#              "twins" = "staging/external_id_maps/dfci-g2c.v1.cleaned.kin0.gz",
#              "out_tsv" = "staging/external_id_maps/dfci-g2c.v1.1KGP_id_map.tsv")

# Load sample metadata, leaving all column names untouched
main.df <- read.table(args$metadata_tsv, header=T, sep="\t",
                      comment.char="", check.names=F)

# Subset to cohort of interest
main.df <- main.df[which(main.df$cohort == args$cohort), 1:2]

# Duplicate external IDs from twins, if optioned
if(!is.null(args$twins)){
  twin.map <- load.kin0.twins(args$twins)
  main.df <- explode.twins(main.df, twin.map)
}

# Write updated manifest to --out-tsv
write.table(main.df, args$out_tsv, col.names=F, row.names=F, quote=F, sep="\t")
