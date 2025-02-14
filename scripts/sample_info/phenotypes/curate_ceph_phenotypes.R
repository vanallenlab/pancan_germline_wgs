#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data downloaded from dbGaP for CEPH


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
G2CR::load.constants("all")

# Declare constants used in variable parsing
sex.map <- c("1" = "male",
             "2" = "female")


##################
# Data Functions #
##################
# Load and populate baseline dataframe from list of included samples
create.df <- function(list.in){
  # Read list of samples
  samples <- sort(unique(read.table(list.in, header=F)[, 1]))

  # Create empty/baseline dataframe with all samples
  data.frame(
    "Sample" = samples,
    "Cohort" = "ceph",
    "reported_sex" = NA,
    "reported_race_or_ethnicity" = "white",
    "age" = NA,
    "birth_year" = NA,
    "vital_status" = NA,
    "age_at_last_contact" = NA,
    "years_to_last_contact" = NA,
    "years_left_censored" = NA,
    "height" = NA,
    "weight" = NA,
    "bmi" = NA,
    "cancer" = "unknown",
    "stage" = NA,
    "metastatic" = NA,
    "grade" = NA,
    "smoking_history" = NA,
    "cancer_icd10" = NA,
    "original_dx" = NA,
    "wgs_tissue" = "blood"
  )
}

# Add sex information from dbGaP
add.sex.info <- function(df, ped.in, sinfo.in){
  # Read & merge dbGaP data
  ped.df <- read.table(ped.in, header=T, blank.lines.skip=TRUE, sep="\t")
  sinfo.df <- read.table(sinfo.in, header=T, blank.lines.skip=TRUE, sep="\t")
  dbgap.df <- merge(ped.df, sinfo.df, sort=F)

  # Fill missing subject IDs
  missing.sid.idxs <- which(dbgap.df$SOURCE_SUBJECT_ID == "")
  dbgap.df$SOURCE_SUBJECT_ID[missing.sid.idxs] <-
    paste("CEPH", dbgap.df$SUBJECT_ID[missing.sid.idxs], sep="_")

  # Map sexes
  sex.v <- remap(dbgap.df$SEX, sex.map)
  names(sex.v) <- dbgap.df$SOURCE_SUBJECT_ID
  df$reported_sex <- sex.v[df$Sample]

  return(df)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate CEPH phenotype data")
parser$add_argument("--sample-ids", metavar=".txt", type="character", required=TRUE,
                    help="List of sample IDs included in analysis")
parser$add_argument("--pedigree-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Pedigree .tsv downloaded from dbGaP")
parser$add_argument("--subject-info-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Subject info .tsv downloaded from dbGaP")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("sample_ids" = "~/scratch/ceph.samples.list",
#              "pedigree_tsv" = "~/Desktop/Collins/VanAllen/UFC/CEPH/CEPH_dbGaP_metadata/phs001872.v1.pht009364.v1.p1.CEPH_Utah_Pedigree.MULTI.txt.gz",
#              "subject_info_tsv" = "~/Desktop/Collins/VanAllen/UFC/CEPH/CEPH_dbGaP_metadata/phs001872.v1.pht009363.v1.p1.CEPH_Utah_Subject.MULTI.txt.gz",
#              "out_tsv" = "~/scratch/ceph.pheno.dev.tsv")

# Generate empty dataframe based on sample IDs
df <- create.df(args$sample_ids)

# Update with sex info from dbGaP
df <- add.sex.info(df, args$pedigree_tsv, args$subject_info_tsv)

# Write to --out-tsv
write.table(df, args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
