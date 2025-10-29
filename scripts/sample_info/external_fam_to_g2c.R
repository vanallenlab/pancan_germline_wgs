#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Map IDs from an input .fam or .ped file for one cohort to G2C study IDs, sexes, and phenotypes


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
parser <- ArgumentParser(description="Convert external .fam/.ped to G2C .fam")
parser$add_argument("--metadata-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Main G2C sample metadata .tsv for all samples")
parser$add_argument("--in-fam", metavar=".fam|.ped", type="character", required=TRUE,
                    help="Original .fam or .ped file for cohort in question")
parser$add_argument("--no-header", action="store_true", default=FALSE,
                    help="Input .fam|.ped file does not have header line")
parser$add_argument("--cohort", metavar="string", type="character", required=TRUE,
                    help="Cohort to evaluate")
parser$add_argument("--out-fam", metavar=".fam", type="character", required=TRUE,
                    help="path to output .fam")
args <- parser$parse_args()

# # DEV:
# args <- list("metadata_tsv" = "~/scratch/dfci-g2c.intake_qc.all.post_qc_batching.non_aou.tsv.gz",
#              "in_fam" = "~/scratch/fam/hgsvc.ped",
#              "no_header" = FALSE,
#              "cohort" = "hgsvc",
#              "out_fam" = "~/scratch/dfci-g2c.hgsvc.fam")
# args <- list("metadata_tsv" = "~/scratch/dfci-g2c.intake_qc.all.post_qc_batching.non_aou.tsv.gz",
#              "in_fam" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/mesa/mesa.fam",
#              "no_header" = FALSE,
#              "cohort" = "mesa",
#              "out_fam" = "~/scratch/dfci-g2c.mesa.fam")

# Load sample metadata, leaving all column names untouched
meta.df <- read.table(args$metadata_tsv, header=T, sep="\t",
                      comment.char="", check.names=F)
meta.df <- meta.df[which(meta.df$cohort == args$cohort), ]
rownames(meta.df) <- meta.df$G2C_id

# Build map of original study IDs to G2C IDs
id.map <- meta.df$G2C_id
names(id.map) <- meta.df$original_id

# Load input .fam|.ped and subset to the first six columns
fam <- G2CR::load.famfile(args$in_fam, header=!args$no_header)

# Prefix all original family designations with cohort for disambiguation
fam$family <- paste(args$cohort, fam$family, sep="_")

# Attempt remapping of all IDs
id.cols <- c("proband", "father", "mother")
fam[, id.cols] <- apply(fam[, id.cols], 2, remap, map=id.map, default.value=0)

# Only retain rows where proband and at least one family member were successfully remapped
fam <- fam[which(apply(fam[, id.cols], 1, function(ids){sum(ids != 0)}) > 1), ]
fam <- fam[which(fam$proband != 0), ]

# Map proband sexes and phenotypes
fam$sex <- remap(meta.df[fam$proband, "batching_sex"],
                 c("male" = 1, "female" = 2), default.value=0)
fam$pheno <- remap(meta.df[fam$proband, "cancer"],
                   c("unknown" = 0, "control" = 1), default.value=2)

# Suffix family IDs with proband G2C ID
# This is necessary to keep a unique mapping of kindred ID vs. trio ID
# (In multiplex families like CEPH this can cause QC issues downstream)
fam$family <- paste(fam$family, fam$proband, sep="_")

# Write updated manifest to --out-fam
write.table(fam, args$out_fam, col.names=F, row.names=F, quote=F, sep="\t")
