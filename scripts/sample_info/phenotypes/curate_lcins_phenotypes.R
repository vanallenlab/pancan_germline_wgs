#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data for LCINS subjects
# Requires supplementary table 1 from Zhang et al., 2021
# https://www.nature.com/articles/s41588-021-00920-0


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
G2CR::load.constants("all")

# Declare constants used in variable parsing
smoking.map <- c("Y" = 1, "N" = 0)
histology.map <- c("other" = "of unspecified histology")
race.map <- c("African" = "black",
              "Asian" = "asian",
              "European" = "white")


##################
# Data Functions #
##################
# Load and clean data from Zhang et al. supp table 1
load.st1 <- function(tsv.in){
  # Read data
  df <- read.table(tsv.in, header=T, sep="\t", quote="")

  # Remap simple columns
  df$Sample <- df$Subject
  df$Cohort <- "lcins"
  df$reported_sex <- tolower(df$gender)
  df$age <- as.numeric(df$age_at_diagnosis)
  df$birth_year <- NA
  df$vital_status <- abs(1-df$death)
  df$age_at_last_contact <- df$age + (df$survival_months / 12)
  df$years_to_last_contact <- df$survival_months / 12
  df$years_left_censored <- 0
  df$bmi <- df$weight <- df$height <- NA
  df$cancer <- "lung"
  df$metastatic <- "unknown"
  df$grade <- df$grade
  df$cancer_icd10 <- NA

  # Per study methods, effectively all samples are caucasian/European
  # Make the blanket assumption that all samples are EUR, which will be
  # updated in a second function based on dbGaP phenotype information (limited).
  # We accept that this will introduce a few (N~4) mislabeled samples, which
  # shouldn't matter because we are going to infer genetic ancestry downstream.
  df$reported_race_or_ethnicity <- "white"

  # Curate cancer diagnosis information
  df$stage <- sub("[A-B]$", "", df$stage)
  df$stage[which(df$stage == "")] <- "unknown"
  df$original_dx <- gsub("[ ]+", "_", sub("s$", "", paste("lung", remap(tolower(df$Histology), histology.map))))

  # Curate smoking data
  df$passive_smoking[which(df$passive_smoking == "")] <- "NA"
  df$smoking_history <- remap(df$passive_smoking, smoking.map)

  # Per the study methods, 48/232 normal WGS samples were sequenced from
  # normal lung tissue. All others were blood. There is no documentation on
  # which 48/232 were from normal lung tissue. The only data available on dbGaP
  # is for the 27 samples from the EAGLE study, which dbGaP lists as lung "body
  # site", but this seems to be in disagreement with prior publications on the
  # EAGLE study, where ~90% of patients have blood normal available.
  #
  # Given that 80% of samples were reported as being from blood, we will label
  # all samples as such. We accept that this will result in mislabeling 20% of
  # the samples. There is a curious subset of 47 (not 48) patients recruited from
  # Nice, France, which is a near match to the 48 normal tissue samples, but
  # absent any confirmation we will not jump to any assumptions
  df$wgs_tissue <- "blood"

  # Return data frame with desired column order
  col.order <- c("Sample", "Cohort", "reported_sex", "reported_race_or_ethnicity",
                 "age", "birth_year", "vital_status", "age_at_last_contact",
                 "years_to_last_contact", "years_left_censored", "height",
                 "weight", "bmi", "cancer", "stage", "metastatic", "grade",
                 "smoking_history", "cancer_icd10", "original_dx", "wgs_tissue")
  return(df[, col.order])
}

# Update sample race/ancestry annotations from dbGaP, where available
add.dbgap.phenos <- function(df, tsv.in){
  # Read dbGaP data
  pheno.df <- read.table(tsv.in, header=T, blank.lines.skip=TRUE, sep="\t")

  # Map onto existing dataframe
  race.v <- remap(pheno.df$RACE, race.map)
  names(race.v) <- pheno.df$SUBJECT_ID
  ovr.idxs <- which(df$Sample %in% names(race.v))
  ovr.ids <- df$Sample[ovr.idxs]
  df$reported_race_or_ethnicity[ovr.idxs] <- race.v[ovr.ids]
  return(df)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate LCINS phenotype data")
parser$add_argument("--supp-table-1-tsv", metavar=".tsv", type="character",
                    required=TRUE, help="Supplementary Table 1 from Zhang et al.")
parser$add_argument("--subject-phenotypes-tsv", metavar=".tsv", type="character",
                    required=TRUE, help="Subject phenotype table from dbGap")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("supp_table_1_tsv" = "~/Desktop/Collins/VanAllen/jackie_younglung/NatGenet_LCINS/NatGenet_LCINS_phenotypes/LCINS.supp_table1.tsv",
#              "subject_phenotypes_tsv" = "~/Desktop/Collins/VanAllen/jackie_younglung/NatGenet_LCINS/NatGenet_LCINS_phenotypes/phs001697.v1.pht010578.v1.p1.c3.EAGLE_Never_Smokers_Subject_Phenotypes.GRU.txt.gz",
#              "out_tsv" = "~/scratch/mesa.pheno.dev.tsv")

# Load and clean data from supplementary table 1
df <- load.st1(args$supp_table_1_tsv)

# Update with reported race info from dbGaP, where available
df <- add.dbgap.phenos(df, args$subject_phenotypes_tsv)

# Write to --out-tsv
write.table(df, args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
