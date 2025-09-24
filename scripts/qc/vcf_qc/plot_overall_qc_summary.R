#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot overall QC summary figures


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
load.constants("all")

# Declare local constants
# TODO: add these


##################
# Data Functions #
##################
# Load summary statistics & organize in sub-dataframes for plotting
load.ss <- function(tsv.in, ref.prefix=NULL){
  ss <- read.table(tsv.in, header=F, sep="\t")
  colnames(ss) <- c("analysis", "measure", "value", "n")

  lapply(c("all", names(var.class.abbrevs)), function(vc){
    count.df <- data.frame()

    ratio.df <- data.frame()

    sb.df <- as.data.frame(rbind(
      ss[which(ss$analysis == gsub("^\\.|\\.$", "",
                                   paste(ref.prefix, "sensitivity.common",
                                         if(vc!="all"){vc}, sep="."))), ],
      ss[which(ss$analysis == gsub("^\\.|\\.$", "",
                                   paste(ref.prefix, "ppv.common",
                                         if(vc!="all"){vc}, sep="."))), ]
    ))

    gb.df <- as.data.frame(rbind(
      ss[which(ss$analysis == paste(vc, "common_hwe", sep=".")), ]
    ))

    inter.df <- data.frame()

    dr.df <- data.frame()

    list("counts" = count.df,
         "ratios" = ratio.df,
         "site_bench" = sb.df,
         "gt_bench" = gb.df,
         "dynamic_ranges" = dr.df)
  })
}


######################
# Plotting Functions #
######################


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot overall QC summary")
parser$add_argument("--stats", metavar=".tsv", type="character",
                    help="File of all QC summary statistics")
parser$add_argument("--site-ref-prefix", metavar="string", type="character",
                    help="String prefix for site-level benchmarking metrics")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV:
# args <- list("stats" = "~/Downloads/dfci-g2c.v1.initial_qc.stats\ 2/dfci-g2c.v1.initial_qc.all_qc_summary_metrics.tsv",
#              "site_ref_prefix" = "gnomad_v4.1",
#              "out_prefix" = "~/scratch/g2c.qc.test")

# TODO: implement this
