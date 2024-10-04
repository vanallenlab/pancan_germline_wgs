#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data downloaded from dbGaP for TOPMed BioMe


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2C, quietly=TRUE)
G2C::load.constants("all")

# Declare constants used in variable parsing
ancestry.map <- c("african-american_african" = "black",
                  "east_south-east_asian" = "east_asian",
                  "european_american" = "white",
                  "hispanic_latin_american" = "latin_american",
                  "multiple_selected" = "other")


##################
# Data Functions #
##################
# Load & reformat BioMe phenotype data
load.phenotypes <- function(tsv.in){
  # Read data
  df <- read.table(tsv.in, header=T, sep="\t", comment.char="#", quote="", na.strings="")

  # Convert simple columns
  df$Sample <- df$SUBJECT_ID
  df$Cohort <- "biome"
  df$reported_sex <- tolower(df$SEX)
  df$age <- as.numeric(apply(df[, grep("age_at_", colnames(df))], 1, min, na.rm=T))
  df$seq_age <- as.integer(df$age_at_dna_blood_draw_wgs)
  df$years_left_censored <- df$seq_age - df$age
  df$birth_year <- NA
  df$age_at_last_contact <- as.numeric(apply(df[, grep("age_at_", colnames(df))], 1, max, na.rm=T))
  anthro.cols <- c("height", "weight", "bmi")
  df[, anthro.cols] <- apply(df[, anthro.cols], 2, as.numeric)
  df$cancer <- "control"
  df$stage <- NA
  df$metastatic <- NA
  df$grade <- NA
  df$cancer_icd10 <- NA
  df$original_dx <- NA
  df$smoking_history <- NA

  # Based on NCBI BioProject information, it appears (?) that all BioMe samples were
  # sequenced from blood-based DNA (e.g., https://www.ncbi.nlm.nih.gov/biosample/13779884)
  # The official BioMe recruitment brochure also indicates only blood samples are collected
  # (https://icahn.mssm.edu/files/ISMMS/Assets/Research/IPM/BioMe-Brochure-English.pdf)
  df$wgs_tissue <- "blood"

  # Infer survival information for individuals whose age at blood draw != oldest reported age
  # Tautologically, these individuals *must* have survived at least this long
  # Treat all other participants as unknown for survival purposes
  df$vital_status <- NA
  followup.idxs <- which(df$age != df$age_at_last_contact)
  df$vital_status[followup.idxs] <- 1
  df$years_to_last_contact <- NA
  df$years_to_last_contact[followup.idxs] <-
    apply(df[followup.idxs, c("age", "age_at_last_contact")], 1, diff)

  # Parse race & ethnicity
  df$ethnicity <- NA
  df$ethnicity[which(!is.na(df$hispanic_subgroup))] <- "hispanic"
  df$race <- remap(tolower(df$ancestry_group), ancestry.map)
  df$roe <- gsub("_NA$", "", apply(df[, c("race", "ethnicity")], 1, paste, collapse="_"))
  df$roe[which(df$roe == "NA" | is.na(df$roe))] <- "unknown"
  df$reported_race_or_ethnicity <- df$roe

  # Return data frame with columns in desired order
  df[, c("Sample", "Cohort", "reported_sex", "reported_race_or_ethnicity", "age",
         "birth_year", "vital_status", "age_at_last_contact", "years_to_last_contact",
         "years_left_censored", "height", "weight", "bmi", "cancer", "stage", "metastatic",
         "grade", "smoking_history", "cancer_icd10", "original_dx", "wgs_tissue")]
}



###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate TOPMed BioMe phenotype data")
parser$add_argument("--in-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="subject_phenotypes.tsv downloaded from dbGaP")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("in_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/topmed/biome/phs001644.v2.pht009946.v2.p2.c1.TOPMed_CCDG_BioME_Subject_Phenotypes.HMB-NPU.txt",
#              "out_tsv" = "~/scratch/biome.pheno.dev.tsv")

# Load and clean data
df <- load.phenotypes(args$in_tsv)

# Write to --out-tsv
write.table(df, args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
