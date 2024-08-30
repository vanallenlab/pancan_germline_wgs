#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data downloaded from dbGaP for GTEx


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2C, quietly=TRUE)
G2C::load.constants("all")

# Declare constants used in variable parsing
sex.map <- c("1" = "male",
             "2" = "female")
race.map <- c("1" = "asian",
              "2" = "black",
              "3" = "white",
              "4" = "native_american")
eth.map <- c("0" = "not_hispanic",
             "1" = "hispanic")


##################
# Data Functions #
##################
# Load & clean subject phenotype information
load.phenotypes <- function(tsv.in){
  # Read data
  df <- read.table(tsv.in, sep="\t", comment.char="#", quote="", header=T,
                   blank.lines.skip=T, check.names=F, fill=T,
                   na.strings=c("", "NA"))
  df <- df[grep("^GTEX-", df$SUBJID), ]

  # Reformat simple columns
  df$Sample <- df$SUBJID
  df$Cohort <- "gtex"
  df$reported_sex <- remap(df$SEX, sex.map, default.value=NA)
  df$age_at_last_contact <- df$age <- as.numeric(df$AGE)
  df$birth_year <- NA
  df$vital_status <- 0
  df$years_to_last_contact <- NA
  df$height <- 2.54 * df$HGHT
  df$weight <- 0.453592 * df$WGHT
  df$bmi <- df$BMI
  df$stage <- NA
  df$metastatic <- NA
  df$grade <- NA
  df$cancer_icd10 <- NA
  df$original_dx <- NA
  df$wgs_tissue <- "unknown"

  # Format race/ethnicity information
  df$race <- remap(df$RACE, race.map, default.value="unknown")
  df$eth <- remap(df$ETHNCTY, eth.map, default.value=NA)
  df$roe <- gsub("_NA$", "", apply(df[, c("race", "eth")], 1, paste, collapse="_"))
  df$reported_race_or_ethnicity <- gsub("^unknown_", "", df$roe)

  # Infer cancer status
  df$cancer <- "control"
  cancer.idxs <- which(df$MHCANCERC == 1 | df$MHCANCER5 == 1)
  df$cancer[cancer.idxs] <- "other"
  df$stage[cancer.idxs] <- "unknown"
  df$metastatic[cancer.idxs] <- "unknown"
  df$grade[cancer.idxs] <- "unknown"

  # Get smoking history
  df$smoke <- as.integer(remap(df$MHSMKSTS, c("Yes" = 1, "No" = 0)))
  df$smokeyears <- as.integer(df$MHSMKYRS)
  df$smokemax <- apply(df[, c("smoke", "smokeyears")], 1, max, na.rm=T)
  df$smokemax[which(is.infinite(df$smokemax))] <- NA
  df$smoking_history <- as.integer(df$smokemax > 0)

  # Return only relevant columns
  keep.cols <- c("Sample", "Cohort", "reported_sex", "reported_race_or_ethnicity",
                 "age", "birth_year", "vital_status", "age_at_last_contact",
                 "years_to_last_contact", "height", "weight", "bmi", "cancer",
                 "stage", "metastatic", "grade", "smoking_history", "cancer_icd10",
                 "original_dx", "wgs_tissue")
  df[, keep.cols]
}

# Add WGS DNA source
add.dna.source <- function(df, tsv.in){
  # Read sample metadata
  s.df <- read.table(tsv.in, header=T, sep="\t", quote="")

  # Infer donor ID
  s.df$Sample <- sapply(s.df$SAMPID, function(x){
    paste(unlist(strsplit(x, split="-"))[1:2], collapse="-")
  })

  # Subset to DNA used for WGS in samples present in df
  s.df <- s.df[which(grepl("WGS", s.df$SMGEBTCHT)
                     & grepl("DNA", s.df$SMMTRLTP)
                     & s.df$Sample %in% df$Sample), ]

  # Deduplicate on donor ID while prioritizing blood over non-blood
  s.df <- s.df[c(which(s.df$SMTS == "Blood"),
                 which(s.df$SMTS != "Blood")), ]
  s.df <- s.df[which(!duplicated(s.df$Sample)), ]
  dna.map <- gsub("[ ]+", "_", tolower(s.df$SMTS))
  names(dna.map) <- s.df$Sample

  # Add DNA source
  df$wgs_tissue <- remap(df$Sample, dna.map)

  return(df)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate GTEx phenotype data")
parser$add_argument("--phenotypes-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Subject_Phenotypes.txt downloaded from dbGaP.")
parser$add_argument("--sample-attributes-tsv", metavar=".tsv", type="character",
                    help="SampleAttributesDS.txt downloaded from GTEx")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("phenotypes_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/gtex/dbgap_phenotypes/phs000424.v9.pht002742.v9.p2.c1.GTEx_Subject_Phenotypes.GRU.txt.gz",
#              "sample_attributes_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/gtex/GTEx_Analysis_2021-02-11_v9_Annotations_SampleAttributesDS.txt",
#              "out_tsv" = "~/scratch/gtex.pheno.dev.tsv")

# Load and clean phenotype data
df <- load.phenotypes(args$phenotypes_tsv)

# Add DNA source
if(!is.null(args$sample_attributes_tsv)){
  df <- add.dna.source(df, args$sample_attributes_tsv)
}

# Write to --out-tsv
col.order <- c("Sample", "Cohort", "reported_sex", "reported_race_or_ethnicity",
               "age", "birth_year", "vital_status", "age_at_last_contact",
               "years_to_last_contact", "height", "weight", "bmi", "cancer",
               "stage", "metastatic", "grade", "smoking_history",
               "cancer_icd10", "original_dx", "wgs_tissue")
write.table(df[, col.order], args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)


