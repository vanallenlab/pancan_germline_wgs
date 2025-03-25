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
cancer.map <- c("breast cancer" = "breast",
                "breast, basal skin" = "breast",
                "colorectal cancer" = "colorectal",
                "melanoma" = "melanoma",
                "prostate cancer" = "prostate",
                "uterine cancer" = "uterus",
                "thyroid" = "thyroid",
                "skin, breast" = "breast",
                "sarcoma on skull" = "sarcoma",
                "colon, melanoma" = "colorectal;melanoma")


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

# Add cancer information from UoU data
add.cancer.info <- function(df, uou.in, sinfo.in, smap.in){
  # Read cancer data
  uou.df <- read.table(uou.in, header=T, sep="\t", quote="")

  # Read sample & subject info and map IDs
  smap <- read.table(smap.in, header=T, blank.lines.skip=TRUE, sep="\t")
  uou.df <- merge(uou.df, smap, by.x="UGRP.Lab.ID", by.y="SAMPLE_ID",
                  all.x=T, all.y=F, sort=F)
  sinfo.df <- read.table(sinfo.in, header=T, blank.lines.skip=TRUE, sep="\t")
  uou.df <- merge(uou.df, sinfo.df, by="dbGaP_Subject_ID", all.x=T, all.y=F,
                  sort=F, suffixes=c("", ".subject"))

  # Map IDs and restrict to samples where unabiguous mapping is possible
  uou.df$Sample <- uou.df$SOURCE_SUBJECT_ID
  uou.df$Sample[which(!uou.df$Sample %in% df$Sample)] <- NA
  sample.ids <- paste("CEPH", uou.df$SUBJECT_ID, sep="_")
  good.sample.ids <- which(sample.ids %in% df$Sample)
  uou.df$Sample[good.sample.ids] <- sample.ids[good.sample.ids]
  uou.df <- uou.df[!is.na(uou.df$Sample), ]

  # Update age info
  age.map <- as.numeric(uou.df$Age.at.last.data.point)
  names(age.map) <- uou.df$Sample
  df$age <- df$age_at_last_contact <- as.numeric(remap(df$Sample, age.map, default.value=NA))
  cancer.age <- sapply(uou.df$SR.cancer.ages, function(v){
    min(as.numeric(gsub(" ", "", unlist(strsplit(v, split=",")))), na.rm=T)
  })
  has.cancer.age <- which(!is.infinite(cancer.age) & !is.na(cancer.age))
  cancer.age.map <- cancer.age[has.cancer.age]
  names(cancer.age.map) <- uou.df$Sample[has.cancer.age]
  df.cancer.age <- as.numeric(remap(df$Sample, cancer.age.map, default.value=NA))
  df$age[which(!is.na(df.cancer.age))] <- df.cancer.age[which(!is.na(df.cancer.age))]
  has.fu <- which(df$age != df$age_at_last_contact)
  df$years_left_censored[has.fu] <- 0
  df$years_to_last_contact[has.fu] <- (df$age_at_last_contact - df$age)[has.fu]
  df$vital_status[has.fu] <- 1

  # Update cancer info
  # Logic applied as follows:
  # 1. Any sample with self-reported or cancer registry Dx is a case
  # 2. Any sample with no recorded cancer data or explicitly "no" known cancer is a control
  # 3. Any sample with disagreement between 1 & 2 is treated as unknown
  #    (This seems to be rare; only two cases with relative-reported cancer but "no" known cancer)
  # 3. All other situations are treated as unknown
  cancer.pos.idx <- grep("CR|RR|SR", uou.df$Cancer.)
  cancer.neg.idx <- union(union(grep("no", uou.df$Cancer.), which(uou.df$Cancer. == "")),
                          grep("no cancer", tolower(uou.df$SR.cancer.type)))
  cancer.unk.idx <- intersect(cancer.pos.idx, cancer.neg.idx)
  uou.df$cancer <- "unknown"
  uou.df$cancer[cancer.pos.idx] <- remap(tolower(uou.df$SR.cancer.type[cancer.pos.idx]),
                                         cancer.map, default.value="other")
  uou.df$cancer[cancer.neg.idx] <- "control"
  uou.df$cancer[cancer.unk.idx] <- "unknown"
  uou.cancer.map <- uou.df$cancer
  names(uou.cancer.map) <- uou.df$Sample
  df$cancer <- remap(df$Sample, uou.cancer.map, default.value = "unknown")

  # Update cancer case-specific metadata
  df.cancer.idx <- which(!df$cancer %in% c("control", "unknown"))
  df[df.cancer.idx, "stage"] <- "unknown"
  df[df.cancer.idx, "metastatic"] <- "unknown"
  df[df.cancer.idx, "grade"] <- "unknown"
  df[df.cancer.idx, "original_dx"] <- gsub(" ", "_", tolower(df[df.cancer.idx, "cancer"]))
  uou.dx.map <- tolower(uou.df$SR.cancer.type[setdiff(cancer.pos.idx, cancer.unk.idx)])
  uou.dx.map <- gsub("_$", "", gsub(" ", "_", gsub(",", "", uou.dx.map)))
  names(uou.dx.map) <- uou.df$Sample[setdiff(cancer.pos.idx, cancer.unk.idx)]
  df$original_dx <- remap(df$Sample, uou.dx.map, default.value=NA)
  df$original_dx[which(df$original_dx == "")] <- NA

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
parser$add_argument("--sample-info-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Sample info .tsv downloaded from dbGaP")
parser$add_argument("--uou-cancer-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Cancer information provided by UoU as .tsv")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("sample_ids" = "~/scratch/ceph.samples.list",
#              "pedigree_tsv" = "~/Desktop/Collins/VanAllen/UFC/CEPH/CEPH_dbGaP_metadata/phs001872.v1.pht009364.v1.p1.CEPH_Utah_Pedigree.MULTI.txt.gz",
#              "subject_info_tsv" = "~/Desktop/Collins/VanAllen/UFC/CEPH/CEPH_dbGaP_metadata/phs001872.v1.pht009363.v1.p1.CEPH_Utah_Subject.MULTI.txt.gz",
#              "sample_info_tsv" = "/Users/ryan/Desktop/Collins/VanAllen/UFC/CEPH/CEPH_dbGaP_metadata/phs001872.v1.pht009365.v1.p1.CEPH_Utah_Sample.MULTI.txt.gz",
#              "uou_cancer_tsv" = "/Users/ryan/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/ceph/CEPH_cancer_data.UoU.jan_2025.tsv",
#              "out_tsv" = "~/scratch/ceph.pheno.dev.tsv")

# Generate empty dataframe based on sample IDs
df <- create.df(args$sample_ids)

# Update with sex info from dbGaP
df <- add.sex.info(df, args$pedigree_tsv, args$subject_info_tsv)

# Update with cancer information from UoU data
df <- add.cancer.info(df, args$uou_cancer_tsv, args$subject_info_tsv,
                      args$sample_info_tsv)

# Write to --out-tsv
write.table(df, args$out_tsv, col.names=T, row.names=F, sep="\t", quote=F)
