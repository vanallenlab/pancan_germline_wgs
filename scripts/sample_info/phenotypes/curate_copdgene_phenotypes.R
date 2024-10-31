#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data downloaded from dbGaP for TOPMed COPDGene


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2C, quietly=TRUE)
G2C::load.constants("all")
sex.map <- c("1" = "male", "2" = "female")
tissue.map <- c("Blood" = "blood", "Nasal epithelial" = "nasal_epithelium")


##################
# Data Functions #
##################
# Load & reformat COPDGene phenotype data
load.phenotypes <- function(sample.tsv, subject.tsv, id.link.tsv){
  # Read sample data
  df <- read.table(sample.tsv, header=T, sep="\t", comment.char="#",
                   quote="", na.strings="")
  df <- df[which(df$ANALYTE_TYPE == "DNA"), ]

  # Read & merge subject data via ID linker tsv
  subj.df <- read.table(subject.tsv, header=T, sep="\t", comment.char="#",
                        quote="", na.strings="")
  link.df <- read.table(id.link.tsv, header=T, sep="\t", comment.char="#",
                        quote="", na.strings="")
  subj.df <- merge(subj.df, link.df, by="dbGaP_Subject_ID", all.x=T, all.y=F, sort=F)
  df <- merge(df, subj.df, by="dbGaP_Sample_ID", all.x=T, all.y=F, sort=F)

  # Convert simple columns
  df$Sample <- df$SAMPLE_ID.x
  df$Cohort <- "copdgene"
  df$reported_sex <- remap(as.character(df$SEX), sex.map)
  df$cancer <- "control"
  df$wgs_tissue <- remap(as.character(df$HISTOLOGICAL_TYPE), tissue.map)
  na.cols <- c("reported_race_or_ethnicity", "age", "birth_year", "vital_status",
               "age_at_last_contact", "years_to_last_contact", "years_left_censored",
               "height", "weight", "bmi", "stage", "metastatic", "grade",
               "smoking_history", "cancer_icd10", "original_dx")
  df[, na.cols] <- NA

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
parser <- ArgumentParser(description="Curate TOPMed COPDGene phenotype data")
parser$add_argument("--sample-attributes", metavar=".tsv", type="character", required=TRUE,
                    help="sample_attributes.tsv downloaded from dbGaP")
parser$add_argument("--topmed-subject-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="TOPMed subject .tsv downloaded from dbGaP")
parser$add_argument("--id-link-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="TOPMed sample-to-subject linker .tsv downloaded from dbGaP")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("sample_attributes" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/copdgene/dbgap_downloads/phs000951.v5.pht005052.v5.p5.c1.TOPMed_WGS_COPDGene_Sample_Attributes.HMB.txt.gz",
#              "topmed_subject_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/copdgene/dbgap_downloads/phs000951.v5.pht005050.v5.p5.TOPMed_WGS_COPDGene_Subject.MULTI.txt.gz",
#              "id_link_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/copdgene/dbgap_downloads/phs000951.v5.pht005051.v5.p5.TOPMed_WGS_COPDGene_Sample.MULTI.txt.gz",
#              "out_tsv" = "~/scratch/copdgene.pheno.dev.tsv")

# Load and clean data
df <- load.phenotypes(args$sample_attributes, args$topmed_subject_tsv, args$id_link_tsv)

# Write to --out-tsv
write.table(df, args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
