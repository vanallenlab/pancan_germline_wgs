#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data for ICGC


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2C, quietly=TRUE)
G2C::load.constants("all")

# Declare constants used in variable parsing
vital.map <- c("alive" = 1,
               "deceased" = 0)
stage.map <- c("STAGE 1" = "I",
               "STAGE 2" = "II", "STAGE 2A" = "II", "STAGE 2B" = "II",
               "STAGE 3" = "III",
               "STAGE 4" = "IV", "STAGE 4S" = "IV",
               "0" = "0",
               "1" = "I", "1A" = "I", "1B" = "I", "1+" = "I",
               "2" = "II", "2A" = "II", "2B" = "II",
               "3" = "III", "3A" = "III", "3B" = "III",
               "4" = "IV", "4A" = "IV", "4B" = "IV",
               "IA" = "I", "IB" = "I", "IAE" = "I",
               "IBE" = "I", "I+" = "I", "IE" = "I", "IC" = "I",
               "IIA" = "II", "IIB" = "II", "IIAE" = "II",
               "IIBE" = "II", "IIC" = "II", "IIE" = "II",
               "IIIA" = "III", "IIIB" = "III", "IIIE" = "III", "IIIC" = "III",
               "IVA" = "IV", "IVB" = "IV", "IVAE" = "IV",
               "METASTATIC STAGE" = "IV",
               "l" = "I", "ll" = "II", "lll" = "III")
pres.map <- c("cryopreservation in dry ice (dead tissue)" = "frozen postmortem",
              "cryopreservation in liquid nitrogen (dead tissue)" = "frozen postmortem",
              "cryopreservation of live cells in liquid nitrogen" = "frozen living",
              "cryopreservation, other" = "frozen",
              "formalin fixed & paraffin embedded" = "FFPE",
              "formalin fixed, unbuffered" = "fixed",
              "other technique" = "")
pres.supp.map <- c("All TCGA cases 1 or 2" = "",
                   "CD3+ tumor cell sorting" = "sorted fresh",
                   "Cryopreservation: liquid nitrogen, dry ice, or other" = "frozen",
                   "fresh into EDTA" = "fresh",
                   "Immersion in RNAlater" = "fresh",
                   "Known to be cryopreservation by method 1, 2, 3, or 4" = "frozen",
                   "Known to be cryopreservation by method 1, 2, or 4" = "frozen")
tissue.supp.map <- c("Non-tumor liver" = "liver",
                     "non tumor liver" = "liver",
                     "Tonsille" = "tonsil",
                     "Saliva" = "saliva",
                     "cerebellum" = "cerebellum",
                     "ENDOMETRIUM" = "endometrium",
                     "stomach" = "stomach",
                     "SKIN" = "skin",
                     "Buccal Cell Normal" = "buccal cells",
                     "spleen" = "spleen",
                     "bone marrow" = "bone marrow",
                     "skin" = "skin",
                     "esophagus" = "esophagus",
                     "bone" = "bone")


##################
# Data Functions #
##################
# Load relevant information from main ICGC donor table
load.donors <- function(tsv.in){
  # Read file
  df <- read.table(tsv.in, header=T, sep="\t", quote="")

  # Convert simple columns
  df$Sample <- df$icgc_donor_id
  df$Cohort <- "icgc"
  df$reported_sex <- df$donor_sex
  df$reported_sex[which(df$reported_sex == "")] <- NA
  df$reported_race_or_ethnicity <- NA
  df$birth_year <- NA
  df$vital_status <- as.numeric(remap(df$donor_vital_status, vital.map,
                                      default.value=NA))
  df$bmi <- df$weight <- df$height <- NA
  df$grade <- "unknown"

  # Curate age information
  age.cols <- grepl("^donor_age_at_", colnames(df))
  df[, age.cols] <- apply(df[, age.cols], 2, function(vals){
    vals <- as.numeric(vals)
    day.idxs <- which(vals > 125)
    vals[day.idxs] <- vals[day.idxs] / 365
    return(vals)
  })
  df$age <- df$donor_age_at_diagnosis
  df$age[which(df$age <= 0 | is.infinite(df$age))] <- NA
  df$fu_age_interval <- (df$donor_age_at_last_followup - df$age) * 365
  fu.t.cols <- grepl("_interval|_time", colnames(df))
  df$max_fu <- apply(df[, fu.t.cols], 1, max, na.rm=T)
  no.fu <- which(df$max_fu <= 0)
  df$years_to_last_contact <- df$max_fu / 365
  df$years_to_last_contact[no.fu] <- NA
  df$years_left_censored <- df$donor_age_at_enrollment - df$donor_age_at_diagnosis
  df$years_left_censored[no.fu] <- NA
  df$years_left_censored[which(is.na(df$years_left_censored) & !is.na(df$years_to_last_contact))] <- 0
  df$age_at_last_contact <- df$age + df$years_to_last_contact
  df$age_at_last_contact[no.fu] <- df$age[no.fu]

  # Curate cancer diagnosis information
  df$cancer_icd10 <- sapply(df$donor_diagnosis_icd10, function(s){
    icd10 <- toupper(as.character(unlist(strsplit(s, split="[, ]"))[1]))
    icd10 <- gsub('^"', "", icd10)
    if(icd10 == "" | icd10 == "N/A" | is.na(icd10)
       | grepl("unknown", icd10, ignore.case=TRUE)){
      NA
    }else{
      icd10 <- gsub("/", ";", icd10)
      if(grepl("^[0-9]", icd10)){paste("C", icd10, sep="")}else{icd10}
    }
  })

  # Curate stage and metastatic status
  stage.cols <- c("donor_tumour_stage_at_diagnosis_supplemental",
                  "donor_tumour_stage_at_diagnosis")
  df[, c("stage", "metastatic")] <- t(apply(df[, stage.cols], 1, function(v){
    # Collate & clean all available info
    for(delim in c("/", ",", "_", " OR ", "-")){
      v <- as.character(unlist(sapply(toupper(v), strsplit, split=delim)))
    }
    v <- sub("( )+$", "", sub("^( )+", "", sub(")", "", sub("(", "", v, fixed=T), fixed=T)))
    v <- remap(v, stage.map)

    # Parse stage info
    # First, check if any obvious stage is reported (after remapping, above)
    met <- stage <- "unknown"
    obvi.stage <- intersect(c("IV", "III", "II", "I", "0"), v)
    if(length(obvi.stage) > 0){
      stage <- obvi.stage[1]
      met <- if(stage == "IV"){1}else{0}
    }else{
      # Otherwise, check if TNM stage is reported, and use T & M to approximate stage
      tnm.v <- gsub("[A-L|O-S|U-Z]", "", v[grep("T[0-4]", v)])
      if(length(tnm.v) > 0){
        # Check if M term is included
        m.v <- tnm.v[grep("M[0-9]", tnm.v)]
        if(length(m.v) > 0){
          met <- sort(sapply(m.v, function(mk){
            parts <- sub("[A-Z]", "", unlist(strsplit(mk, split="M")))
            as.numeric(parts[length(parts)])
          }), decreasing=TRUE)[1]
          if(met == 1){stage <- "IV"}
        }
        # Check if T term is included
        m.t <- tnm.v[grep("T[0-9]", tnm.v)]
        if(length(m.t) > 0){
          stage <- remap(sort(sapply(m.t, function(tk){
            parts <- as.numeric(unlist(strsplit(tk, split="[A-Z]")))
            parts[which(!is.na(parts))][1]
          }), decreasing=TRUE)[1], stage.map, default.value="unknown")
          if(stage == "IV"){met <- 1}
        }
      }
    }
    c(stage, met)
  }))

  return(df)
}

# Annotate cancer diagnoses and G2C cancer types per sample
# This requires a three-column .tsv input with ICGC project code,
# description, and G2C cancer type. This was sourced from:
# http://www.innovebioinfo.com/Database/SMDB/sample_file/ICGC.html
# https://github.com/SimonHensel/health-data-lake?tab=readme-ov-file
annotate.cancers <- function(df, tsv.in){
  # Read cancer map
  c.map <- read.table(tsv.in, header=T, sep="\t")
  colnames(c.map)[2] <- "original_dx"
  c.map$original_dx <- gsub("[ ]+", "_", tolower(c.map$original_dx))

  # Left outer join to map cancers to project IDs
  merge(df, c.map, all.x=T, all.y=F, sort=F)
}

# Annotate smoking history
annotate.smoking <- function(df, tsv.in){
  # Read exposure data
  ex.df <- read.table(tsv.in, header=T, sep="\t", quote="")

  # Build smoking map
  smoke.cols <- c("tobacco_smoking_intensity", "tobacco_smoking_history_indicator")
  smoke.map <- apply(ex.df[, smoke.cols], 1, function(v){
    si <- if(is.na(v[1])){-999}else{as.numeric(v[1])}
    sh <- as.character(v[2])
    if(si > 0){
      return(1)
    }else if(sh == "Lifelong non-smoker (<100 cigarettes smoked in lifetime)"){
      return(0)
    }else if(sh != "Smoking history not documented" & !is.na(sh) & sh != ""){
      return(1)
    }else{
      return(NA)
    }
  })
  names(smoke.map) <- ex.df$icgc_donor_id

  # Map smoking history
  df$smoking_history <- remap(df$Sample, smoke.map, default.value=NA)
  return(df)
}

# Annotate DNA source tissue
annotate.tissue <- function(df, tsv.in){
  # Read specimen.tsv and subset to normal samples
  spec <- read.table(tsv.in, header=T, sep="\t")
  norm.idxs <- apply(spec[, c("specimen_type", "specimen_type_other")], 1,
                     function(v){
                       any(grepl("normal", v, ignore.case=T)) &
                         !any(grepl("tumour", v, ignore.case=T))
                     })
  spec <- spec[norm.idxs, ]

  # Reformat specimen information
  spec$tissue <- tolower(gsub(" derived", "", gsub("^Normal - ", "", spec$specimen_type)))
  no.tissue.idx <- which(spec$tissue %in% c("ebv immortalized", "other", "solid tissue",
                                            "tissue adjacent to primary"))
  spec$tissue_backup <- NA
  spec$tissue_backup[no.tissue.idx] <-
    remap(spec$specimen_type_other[no.tissue.idx], tissue.supp.map, default.value=NA)
  spec$final_tissue <- apply(spec[, c("tissue_backup", "tissue")], 1, function(v){
    head(v[which(!is.na(v))], 1)
  })
  spec$storage <- remap(spec$specimen_processing, pres.map)
  no.storage.idx <- which(spec$storage == "")
  spec$storage[no.storage.idx] <- remap(spec$specimen_processing_other[no.storage.idx],
                                        pres.supp.map)
  tissue.map <- apply(spec[, c("storage", "final_tissue")], 1, paste, collapse=" ")
  tissue.map <- gsub("^_", "", tolower(gsub(" ", "_", tissue.map)))
  names(tissue.map) <- spec$icgc_donor_id

  # Map specimen info onto dataframe
  df$wgs_tissue <- remap(df$Sample, tissue.map, default.value="unknown")
  df$wgs_tissue[which(is.na(df$wgs_tissue))] <- "unknown"
  return(df)
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate ICGC phenotype data")
parser$add_argument("--donors-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="donor.tsv.gz downloaded from ICGC")
parser$add_argument("--project-cancer-map", metavar=".tsv", type="character",
                    help="three-column .tsv mapping ICGC project code to G2C cancer",
                    required=TRUE)
parser$add_argument("--specimen-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="specimen.tsv.gz downloaded from ICGC")
parser$add_argument("--exposure-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="donor_exposure.tsv.gz downloaded from ICGC")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("donors_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/icgc/icgc_wgs_download_may26_2023/icgc_donor_metadata_release_28_may_2023/donor.all_projects.tsv.gz",
#              "project_cancer_map" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/icgc/icgc_project_to_g2c_map.tsv",
#              "specimen_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/icgc/icgc_wgs_download_may26_2023/icgc_donor_metadata_release_28_may_2023/specimen.all_projects.tsv.gz",
#              "exposure_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/icgc/icgc_wgs_download_may26_2023/icgc_donor_metadata_release_28_may_2023/donor_exposure.all_projects.tsv.gz",
#              "out_tsv" = "~/scratch/icgc.pheno.dev.tsv")

# Load main donor table
df <- load.donors(args$donors_tsv)

# Annotate with cancer diagnosis and G2C cancer
df <- annotate.cancers(df, args$project_cancer_map)

# Add smoking history
df <- annotate.smoking(df, args$exposure_tsv)

# Annotate DNA tissue
df <- annotate.tissue(df, args$specimen_tsv)

# Write to --out-tsv
col.order <- c("Sample", "Cohort", "reported_sex", "reported_race_or_ethnicity",
               "age", "birth_year", "vital_status", "age_at_last_contact",
               "years_to_last_contact", "years_left_censored", "height",
               "weight", "bmi", "cancer", "stage", "metastatic", "grade",
               "smoking_history", "cancer_icd10", "original_dx", "wgs_tissue")
write.table(df[, col.order], args$out_tsv, col.names=T,
            row.names=F, sep="\t", quote=F)
