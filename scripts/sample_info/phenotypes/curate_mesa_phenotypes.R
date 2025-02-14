#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data downloaded from dbGaP for TOPMed MESA


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
G2CR::load.constants("all")

# Declare constants used in variable parsing
sex.map <- c("0" = "female",
             "1" = "male")
race.map <- c("1" = "white",
              "2" = "east_asian",
              "3" = "black",
              "4" = "hispanic")
cancer.map <- c("0" = "control",
                "1" = "other",
                "9" = "control")


##################
# Data Functions #
##################
# Load and clean data from initial study intake exam
load.intake.exam <- function(tsv.in, col.suffix="1"){
  # Read data
  df <- read.table(tsv.in, sep="\t", comment.char="#", header=T, quote="", na.strings="")

  # Convert simple columns
  df$Sample <- paste("ssi", df$sidno, sep="_")
  df$Cohort <- "mesa"
  df$reported_sex <- remap(df[, paste("gender", col.suffix, sep="")], sex.map)
  df$reported_race_or_ethnicity <- remap(df[, paste("race", col.suffix, "c", sep="")],
                                         race.map, default.value="unknown")
  df$vital_status <- NA
  df$age_at_last_contact <- df$age <- as.numeric(df[, paste("age", col.suffix, "c", sep="")])
  df$birth_year <- NA
  df$years_left_censored <- df$years_to_last_contact <- NA
  df$height <- as.numeric(df[, paste("htcm", col.suffix, sep="")])
  df$weight <- 0.453592 * as.numeric(df[, paste("wtlb", col.suffix, sep="")])
  df$bmi <- df$weight / ((df$height / 100)^2)
  df$stage <- NA
  df$metastatic <- NA
  df$grade <- NA
  df$smoking_history <- as.integer(as.integer(df[, paste("cig", col.suffix, "c", sep="")]) > 0)
  df$cancer_icd10 <- NA
  df$original_dx <- NA
  df$wgs_tissue <- "blood"

  # Get cancer labels
  df$prostate <- remap(df[, paste("prostcn", col.suffix, sep="")],
                       c("1" = "prostate"), default.value=NA)
  df$breast <- remap(df[, paste("brstcn", col.suffix, sep="")],
                     c("1" = "breast"), default.value=NA)
  df$lung <- remap(df[, paste("lungcn", col.suffix, sep="")],
                   c("1" = "lung"), default.value=NA)
  df$crc <- remap(df[, paste("coloncn", col.suffix, sep="")],
                  c("1" = "colorectal"), default.value=NA)
  df$cstr <- apply(df[, c("prostate", "breast", "lung", "crc")], 1,
                   function(v){paste(v[which(!is.na(v))], collapse=";")})
  cstr.hits <- which(df$cstr != "")
  df$cancer <- remap(df[, paste("cancer", col.suffix, sep="")], cancer.map)
  df$cancer[cstr.hits] <- df$cstr[cstr.hits]
  any.cancer <- which(df$cancer != "control" & !is.na(df$cancer))
  df$stage[any.cancer] <- "unknown"
  df$metastatic[any.cancer] <- "unknown"
  df$grade[any.cancer] <- "unknown"
  na.cancer <- is.na(df$cancer)
  df$cancer[na.cancer] <- "unknown"

  # Only return relevant columns
  keep.cols <- c("Sample", "Cohort", "reported_sex", "reported_race_or_ethnicity",
                 "age", "birth_year", "vital_status", "age_at_last_contact",
                 "years_to_last_contact", "years_left_censored", "height",
                 "weight", "bmi", "cancer", "stage", "metastatic", "grade",
                 "smoking_history", "cancer_icd10", "original_dx", "wgs_tissue")
  df[, keep.cols]
}

# Update phenotype dataframe with longitudinal follow-up, where available
add.followup <- function(df, tsv.in){
  # Read data
  fu.df <- read.table(tsv.in, sep="\t", comment.char="#", header=T, quote="", na.strings="")

  # Restrict to samples present in df
  fu.df$Sample <- paste("ssi", fu.df$sidno, sep="_")
  fu.df <- fu.df[which(fu.df$Sample %in% df$Sample), ]
  rownames(fu.df) <- fu.df$Sample

  # Infer last known age
  age.v <- fu.df[, grep("age[2-9]c", colnames(fu.df))]
  names(age.v) <- rownames(fu.df)
  age.v <- age.v[which(!is.na(age.v))]
  df$age_at_last_contact <- as.numeric(remap(df$Sample, age.v))
  has.fu <- which(!is.na(df$age_at_last_contact))
  df$vital_status[has.fu] <- 1
  df$years_left_censored[has.fu] <- 0

  # Update inferred survival time
  df$years_to_last_contact[has.fu] <- (df$age_at_last_contact - df$age)[has.fu]

  return(df)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate TOPMed MESA phenotype data")
parser$add_argument("--first-exam-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Earliest MESA main exam downloaded from dbGaP")
parser$add_argument("--followup-tsv", metavar=".tsv", type="character", action="append",
                    help=paste("Follow-up MESA main exam downloaded from dbGaP.",
                               "If multiple exams are provided, they must be",
                               "specified in chronological order."))
parser$add_argument("--family-exam-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Main exam for MESA family participants downloaded from dbGaP")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("first_exam_tsv" = "~/Desktop/Collins/VanAllen/RAS_projects/RAS_WGS/dbGaP/MESA_phs001416.v2.p1/manifests/redownload_dec19_2023/MESA_tsvs/phs000209.v13.pht001116.v10.p3.c1.MESA_Exam1Main.HMB.txt.gz",
#              "followup_tsv" = c("~/Desktop/Collins/VanAllen/RAS_projects/RAS_WGS/dbGaP/MESA_phs001416.v2.p1/manifests/redownload_dec19_2023/MESA_tsvs/phs000209.v13.pht001118.v8.p3.c1.MESA_Exam2Main.HMB.txt.gz",
#                                  "~/Desktop/Collins/VanAllen/RAS_projects/RAS_WGS/dbGaP/MESA_phs001416.v2.p1/manifests/redownload_dec19_2023/MESA_tsvs/phs000209.v13.pht001119.v8.p3.c1.MESA_Exam3Main.HMB.txt.gz",
#                                  "~/Desktop/Collins/VanAllen/RAS_projects/RAS_WGS/dbGaP/MESA_phs001416.v2.p1/manifests/redownload_dec19_2023/MESA_tsvs/phs000209.v13.pht001120.v10.p3.c1.MESA_Exam4Main.HMB.txt.gz",
#                                  "~/Desktop/Collins/VanAllen/RAS_projects/RAS_WGS/dbGaP/MESA_phs001416.v2.p1/manifests/redownload_dec19_2023/MESA_tsvs/phs000209.v13.pht003091.v3.p3.c1.MESA_Exam5Main.HMB.txt.gz"),
#              "family_exam_tsv" = "~//Desktop/Collins/VanAllen/RAS_projects/RAS_WGS/dbGaP/MESA_phs001416.v2.p1/manifests/PhenotypeFiles/phs000209.v13.phenotypes/phs000209.v13.pht001121.v3.p3.c1.MESA_FamilyExamMain.HMB.txt.gz",
#              "out_tsv" = "~/scratch/mesa.pheno.dev.tsv")

# Load and clean data at intake
df <- load.intake.exam(args$first_exam_tsv)
fam.df <- load.intake.exam(args$family_exam_tsv, col.suffix="f")
df <- as.data.frame(rbind(df, fam.df[which(!(fam.df$Sample %in% df$Sample)), ]))

# Add longitudinal follow-up information
if(!is.null(args$followup_tsv)){
  for(followup.tsv in args$followup_tsv){
    df <- add.followup(df, followup.tsv)
  }
}

# Write to --out-tsv
col.order <- c("Sample", "Cohort", "reported_sex", "reported_race_or_ethnicity",
               "age", "birth_year", "vital_status", "age_at_last_contact",
               "years_to_last_contact", "years_left_censored", "height",
               "weight", "bmi", "cancer", "stage", "metastatic", "grade",
               "smoking_history", "cancer_icd10", "original_dx", "wgs_tissue")
write.table(df[, col.order], args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
