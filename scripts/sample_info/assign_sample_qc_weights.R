#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Assign sample selection probabilities for QC


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
# Parse command line arguments and opions
parser <- ArgumentParser(description="Compute sampling probabilities")
parser$add_argument("--qc-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Intake QC .tsv")
parser$add_argument("--pass-column", metavar="string", type="character",
                    help=paste("Intake QC .tsv will be filtered according to",
                               "boolean (or boolean-coercible) 'true' in",
                               "this column prior to plotting. Column values",
                               "are not case sensitive. Can be specified",
                               "multiple times"), action="append")
parser$add_argument("--out-tsv", metavar="path", type="character", required=TRUE,
                    help="Path to output .tsv of sample IDs and sampling probabilities")
args <- parser$parse_args()

# # DEV:
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.post_qc_batching.tsv.gz",
#              "pass_column" = c("global_qc_pass", "batch_qc_pass"),
#              "out_tsv" = "~/scratch/dfci-g2c.qc_probabilities.test.tsv")

# Load data
qc.df <- load.sample.qc.df(args$qc_tsv)

# Handle --pass-column behavior, if optioned
if(!is.null(args$pass_column) & length(args$pass_column) > 0){
  # Collect list of row indexes failing any of --pass-column
  fail.rows <- c()
  for(c.name in args$pass_column){
    if(c.name %in% colnames(qc.df)){
      if(is.numeric(qc.df[, c.name])){
        fail.rows <- c(fail.rows, which(!as.logical(qc.df[, c.name])))
      }else{
        fail.rows <- c(fail.rows, which(!as.logical(toupper(as.character(qc.df[, c.name])))))
      }
    }else{
      stop(paste("Column name '", c.name, "', specified as --pass-column,",
                 "but this column could not be located in --qc-tsv. Exiting."))
    }
  }
  fail.rows <- sort(unique(fail.rows))
  if(length(fail.rows) > 0){
    qc.df <- qc.df[-fail.rows, ]
  }
}

# Assign each sample to a stratum
strata <- apply(qc.df[, c("cohort", "batching_sex", "intake_qc_pop")], 1, paste, collapse="_")

# Compute sample weights to ensure representative sampling across cohort, sex, and ancestry
strata.probs <- 1 / table(strata)
sample.probs <- as.numeric(remap(as.character(strata), strata.probs))
sample.weights <- (nrow(qc.df) / sum(sample.probs)) * sample.probs
names(sample.weights) <- names(strata)

# Write out sampling weights to --out-tsv
out.df <- data.frame("#G2C_ID" = names(sample.weights),
                     "weight" = as.numeric(sample.weights),
                     check.names=F)
write.table(out.df, args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
