#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data for DFCI-Proactive


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(readxl, quietly=TRUE)
require(G2CR, quietly=TRUE)
G2CR::load.constants("all")

# Declare constants used in variable parsing
cohort.map <- c("TRUE" = "proactive-core", "FALSE" = "proactive-other")
race.map <- c("african american" = "black",
              "american indian or alaska native" = "native_american",
              "arab" = "middle_eastern",
              "black or african american" = "black",
              "chinese" = "east_asian",
              "declined" = NA,
              "n/a" = NA,
              "na" = NA,
              "native hawaiian or other islander" = "native_american",
              "not reported" = NA,
              "not specified" = NA,
              "some other race" = "other",
              "unknown" = NA,
              "unknown - declined" = NA)
eth.map <- c("h" = "hispanic",
             "hispanic" = "hispanic",
             "hispanic or latino" = "hispanic",
             "hispanic/latino" = "hispanic",
             "n" = "not_hispanic",
             "non - hispanic" = "not_hispanic",
             "non hispanic" = "not_hispanic",
             "non-hispanic" = "not_hispanic",
             "not hispanic" = "not_hispanic",
             "other" = "other",
             "y" = "hispanic",
             "yes" = "hispanic")
cancer.map <- c("Bladder" = "bladder",
                "Breast" = "breast",
                "Cholangiocarcinoma" = "other",
                "CPOP Referral" = "other",
                "Endometrium" = "uterus",
                "Gastric" = "stomach",
                "Genetics Referral" = "other",
                "Glioblastoma" = "cns",
                "Lung" = "lung",
                "Melanoma" = "melanoma",
                "Mesothelioma" = "lung",
                "Multiple Myeloma" = "other",
                "NUT Carcinoma" = "other",
                "Other" = "other",
                "Pediatric" = "other",
                "Renal" = "kidney",
                "Sarcoma" = "sarcoma",
                "Therapy Associated Polyposis (TAP)" = "other",
                "Thyroid" = "thyroid",
                "Young Colorectal Cancer (CRC)" = "colorectal",
                "Young Lung (dx < 46yrs)" = "lung",
                "Young Tongue" = "oral_cavity")
yl.sex.map <- c("M" = "male", "F" = "female")
yl.race.map <- c("african" = "black",
                 "brazilian" = "american",
                 "caucasian" = "white",
                 "pacific_islander" = "american")
yl.hist.map <- c("lcne" = "large_cell_neuroendocrine_tumor",
                 "luad" = "adenocarcinoma",
                 "neuroend" = "neuroendocrine_tumor",
                 "scc" = "squamous_cell_carcinoma",
                 "sclc" = "small_cell_cancer",
                 "unknown" = "cancer")
yl.stage.map <- c("1" = "I",
                  "2" = "II",
                  "3" = "III",
                  "4" = "IV")


##################
# Data Functions #
##################
# Load and curate a single --meta-xlsx
load.main.meta <- function(tsv.in){
  # Read excel file
  meta <- as.data.frame(readxl::read_xlsx(tsv.in))

  # Handle duplicate columns
  col.names <- gsub("...[0-9]+$", "", colnames(meta))
  keep.idxs <- which(!duplicated(col.names))
  meta <- meta[, keep.idxs]
  colnames(meta) <- col.names[keep.idxs]

  # Rename columns used for collapsing
  meta$sample_id <- meta$`Broad Specimen ID for WGS and RNA (Terra: sample_id)`
  meta$batch <- meta$`Batch number in which specimen was sent to the Broad for WGS/RNA (Terra: sample_set_id)`

  return(meta)
}

# Read, merge, and deduplicate all --meta-xlsx
load.main.metas <- function(tsvs.in){
  # Read first df
  meta <- load.main.meta(tsvs.in[1])

  # If more than one --meta-xlsx are provided, iterative append rows for
  # non-duplicate samples from each
  if(length(tsvs.in) > 1){
    for(new.tsv in tsvs.in[-1]){
      meta.new <- load.main.meta(new.tsv)

      # Identify non-duplicate samples based on MRN + Broad sample ID
      keep.ridx <- which(!(meta.new$`OncDRS DFCI MRN` %in% meta$`OncDRS DFCI MRN`)
                         & !(meta.new$sample_id %in% meta$sample_id)
                         & !is.na(meta.new$sample_id))

      # If any new samples are found, merge into main manifest
      if(length(keep.ridx) > 0){
        meta.new <- meta.new[keep.ridx, ]

        # Unify column names
        all.cols <- union(colnames(meta), colnames(meta.new))
        meta[, setdiff(all.cols, colnames(meta))] <- NA
        meta.new[, setdiff(all.cols, colnames(meta.new))] <- NA

        meta <- rbind(meta, meta.new)
      }
    }
  }

  return(meta)
}

# Add missing samples from Terra metadata
add.terra.info <- function(df, p.in, s.in){
  # Read & merge Terra participant.tsv and sample.tsv
  p.df <- read.table(p.in, header=T, sep="\t", comment.char="")
  s.df <- read.table(s.in, header=T, sep="\t", comment.char="")
  terra.df <- merge(p.df, s.df, by="collaborator_participant_id", all=T, sort=F)

  # Reformat to contain minimal necessary columns
  terra.df$`OncDRS Sex` <- terra.df$gender
  terra.df$sample_id <- gsub("-N-[0-9]+$", "", gsub("-DNA[1-9]$", "", terra.df$sample))
  missing.sid <- which(terra.df$sample_id == "" | is.na(terra.df$sample_id))
  terra.df$sample_id[missing.sid] <- terra.df$collaborator_participant_id[missing.sid]
  terra.df[, setdiff(colnames(terra.df), colnames(df))] <- NULL

  # Fill NAs for columns from main metadata not present in Terra manifest
  # (Note: this will be almost everything except for sex)
  terra.df[, setdiff(colnames(df), colnames(terra.df))] <- NA

  # Identify samples only present in Terra participants and append to main df
  # Also reorders columns to match main manifest df
  add.df <- terra.df[which(!terra.df$sample_id %in% df$sample_id), colnames(df)]

  # Add new samples and return combined dataframe
  as.data.frame(rbind(df, add.df))
}

# Add missing samples based on a list of known processed sample IDs
add.missing.samples <- function(df, tsv.in){
  # Read list of missing smaples
  all.ids <- sort(unique(read.table(tsv.in, header=F)[, 1]))

  # Find sample IDs not already in main metadata
  ids.to.add <- setdiff(all.ids, df$sample_id)

  # Do nothing if no missing IDs are found
  if(length(ids.to.add) == 0){
    return(df)
  }

  # Create dummy dataframe for missing samples
  add.df <- data.frame("sample_id" = ids.to.add)
  add.df[, setdiff(colnames(df), "sample_id")] <- NA

  # Add new samples and return combined dataframe
  as.data.frame(rbind(df, add.df))
}

# Curate values in merged main metadata from load.main.metas
curate.main.meta <- function(df){
  # Curate simple columns
  df$proca_cohort <- df$Cohort
  df$Cohort <- remap(grepl("^PROCA", df$sample_id), cohort.map)
  df$Sample <- gsub("\\.$", "", df$sample_id)
  df$reported_sex <- remap(tolower(df$`OncDRS Sex`), c("unknown" = NA))
  df$age <- as.numeric(df$`Age at Consent`)
  df$birth_year <- sapply(df$`OncoDRS Date of Birth`, function(s){
    as.numeric(unlist(strsplit(as.character(s), split="-", fixed=T))[1])
  })
  df$wgs_tissue <- "blood"

  # Fill columns for which we don't have any data
  na.cols <- c("vital_status", "age_at_last_contact", "years_to_last_contact",
               "years_left_censored", "height", "weight", "bmi",
               "smoking_history", "cancer_icd10")
  df[, na.cols] <- NA
  unk.cols <- c("stage", "metastatic", "grade")
  df[, unk.cols] <- "unknown"

  # Curate race/ethnicity data
  df$race <- remap(tolower(df$`OncDRS Race`), race.map)
  df$eth <- remap(tolower(df$`OncDRS Hispanic Ethnicity per NAACCR`),
                  eth.map, default.value=NA)
  df$reported_race_or_ethnicity <- apply(df[, c("race", "eth")], 1, function(v){
    v <- as.character(v[which(!is.na(v))])
    if(length(v) == 0){NA}else{tolower(paste(v, collapse="_"))}
  })

  # Curate cancer diagnosis information
  df$cancer <- remap(df$proca_cohort, cancer.map)
  df$cancer[grep("^[E|G]OYLC-", df$Sample)] <- "lung"
  df$cancer[grep("^RC_", df$Sample)] <- "melanoma"
  df$cancer[grep("_GL_", df$Sample)] <- "other"
  df$cancer[which(is.na(df$cancer))] <- "other"
  df$original_dx <- apply(df[, grep("Diagnosis", colnames(df))], 1, function(v){
    v <- as.character(v[which(!is.na(v))])
    if(length(v) == 0){NA}else{gsub("'", "", gsub("[ ]+", "_", tolower(paste(v, collapse="_"))))}
  })

  return(df)
}

# Update metadata for young lung cases from YL
update.yl.meta <- function(df, tsv.in){
  # Read young lung metadata
  yl.df <- read.table(tsv.in, header=T, sep="\t", comment.char="")

  # Drop PC, PRS, pollution columns (not needed)
  drop.col.idxs <- unique(c(grep("PC", colnames(yl.df)),
                            grep("score$", colnames(yl.df)),
                            grep("[0-9]yr$", colnames(yl.df))))
  yl.df <- yl.df[, -drop.col.idxs]

  # Subset to samples present in main metadata and reorder to match
  main.yl.idx <- which(df$sample_id %in% yl.df$sample_id)
  rownames(yl.df) <- yl.df$sample_id
  yl.df <- yl.df[df$sample_id[main.yl.idx], ]

  # Update demographic & descriptive metadata
  df$reported_sex[main.yl.idx] <- remap(yl.df$sex, yl.sex.map)
  df$reported_race_or_ethnicity[main.yl.idx] <-
    remap(gsub("[ ]+$", "", tolower(yl.df$race_reported)), yl.race.map)
  df$height[main.yl.idx] <- as.numeric(yl.df$Height_cm)
  df$weight[main.yl.idx] <- as.numeric(yl.df$Weight_kg)
  df$bmi[main.yl.idx] <- as.numeric(yl.df$BMI)

  # Update date/time-based metadata
  df$age[main.yl.idx] <- as.numeric(yl.df$age_dx)
  yl.df$DOB[which(yl.df$DOB == "")] <- NA
  yl.df$yob <- as.numeric(format(as.Date(as.character(yl.df$DOB), format="%m/%d/%y", origin = "1904-01-01"), "%Y"))
  impossible.yob <- which(yl.df$yob > 2023)
  yl.df$yob[impossible.yob] <- yl.df$yob[impossible.yob] - 100
  df$birth_year[main.yl.idx] <- yl.df$yob
  df$vital_status[main.yl.idx] <- as.numeric(yl.df$Survival_status)
  yl.df$age_at_last_contact <- as.numeric(yl.df$age_dx) + as.numeric(yl.df$Survival_years)
  yl.df$age_at_last_contact[is.na(yl.df$age_at_last_contact)] <- yl.df$age_dx[is.na(yl.df$age_at_last_contact)]
  df$age_at_last_contact[main.yl.idx] <- yl.df$age_at_last_contact
  df$years_to_last_contact[main.yl.idx] <- as.numeric(yl.df$Survival_years)
  df$years_left_censored[main.yl.idx] <- 0

  # Update cancer diagnostic metadata
  df$cancer[main.yl.idx] <- "lung"
  yl.df$stage <- as.numeric(yl.df$stage_dx_whole)
  yl.df$stage[is.na(yl.df$stage)] <- "unknown"
  df$stage[main.yl.idx] <- remap(as.character(yl.df$stage), yl.stage.map)
  yl.df$met <- as.numeric(yl.df$stage == 4)
  yl.df$met[which(yl.df$stage == "unknown")] <- "unknown"
  df$metastatic[main.yl.idx] <- yl.df$met
  df$smoking_history[main.yl.idx] <- as.numeric(yl.df$smok_status > 0)
  yl.df$orig_dx <- paste("lung", remap(tolower(yl.df$histology), yl.hist.map), sep="_")
  df$original_dx[main.yl.idx] <- yl.df$orig_dx

  # Ensure no missing tissue source
  df$wgs_tissue[which(is.na(df$wgs_tissue))] <- "unknown"

  # Return updated data
  return(df)
}

# Deduplicate phenotype data.frame based on fewest number of missing variables
dedup.df <- function(df){
  n.missing <- apply(df, 1, function(v){length(which(is.na(v)))})
  df <- df[order(n.missing), ]
  df <- df[!duplicated(df$Sample), ]
  df[order(df$Sample), ]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate DFCI-Proactive phenotype data")
parser$add_argument("--meta-xlsx", metavar=".xlsx", type="character",
                    required=TRUE, action="append",
                    help=paste(".xlsx of Proactive sample information. Must be",
                               "specified once and can be supplied multiple times.",
                               "Should be provided in reverse chronological order",
                               "(newest first)."))
parser$add_argument("--participant-tsv", metavar=".tsv", type="character",
                    required=TRUE, help="participant.tsv from Terra")
parser$add_argument("--sample-tsv", metavar=".tsv", type="character",
                    required=TRUE, help="sample.tsv from Terra")
parser$add_argument("--processed-samples", metavar=".txt", type="character",
                    required=TRUE, help="list of all sample IDs processed to date")
parser$add_argument("--young-lung-meta", metavar=".tsv", type="character", required=TRUE,
                    help="clinical metadata for young lung cases provided by JL")
parser$add_argument("--out-dir", metavar="path", type="character", required=TRUE,
                    help="path to output directory")
args <- parser$parse_args()

# # DEV:
# args <- list("meta_xlsx" = c("~/Partners HealthCare Dropbox/Ryan Collins/PROCA_DATA/PROCA IDs for PROACTIVE Data_15Mar24.xlsx",
#                              "~/Partners HealthCare Dropbox/Ryan Collins/PROCA_DATA/TOP PROACTIVE Cohort_PROCA IDs_4June24.xlsx",
#                              "~/Partners HealthCare Dropbox/Ryan Collins/PROCA_DATA/PROCA IDs for PROACTIVE Data_8Feb24.xlsx",
#                              "~/Partners HealthCare Dropbox/Ryan Collins/PROCA_DATA/PROCA IDs for PROACTIVE Data_7Nov23.xlsx",
#                              "~/Partners HealthCare Dropbox/Ryan Collins/PROCA_DATA/PROCA IDs for PROACTIVE Data_21Nov23.xlsx",
#                              "~/Partners HealthCare Dropbox/Ryan Collins/PROCA_DATA/PROACTIVE All Consented 15Aug23_for Jackie.xlsx"),
#              "participant_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/proactive/GARBER-PROACTIVE_WGS_terra_metadata_sept24_2024/participant.tsv",
#              "sample_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/proactive/GARBER-PROACTIVE_WGS_terra_metadata_sept24_2024/sample.tsv",
#              "processed_samples" = "~/scratch/proactive.all.samples.list",
#              "young_lung_meta" = "~/Desktop/Collins/VanAllen/jackie_younglung/younglung_metadata/combined_data_ALL.tsv",
#              "out_dir" = "~/scratch")

# Read, merge, and deduplicate all --meta-xlsx provided
df <- load.main.metas(args$meta_xlsx)

# Add missing samples from Terra metadata
df <- add.terra.info(df, args$participant_tsv, args$sample_tsv)

# Add missing samples from --samples-list
df <- add.missing.samples(df, args$processed_samples)

# Curate main metadata
df <- curate.main.meta(df)

# Update missing metadata from young lung manifest
df <- update.yl.meta(df, args$young_lung_meta)

# Deduplicate by sample ID based on least number of missing variables
df <- dedup.df(df)

# Write to --out-tsv separately for proactive-core and proactive-other
col.order <- c("Sample", "Cohort", "reported_sex", "reported_race_or_ethnicity",
               "age", "birth_year", "vital_status", "age_at_last_contact",
               "years_to_last_contact", "years_left_censored", "height",
               "weight", "bmi", "cancer", "stage", "metastatic", "grade",
               "smoking_history", "cancer_icd10", "original_dx", "wgs_tissue")
write.table(df[which(df$Cohort == "proactive-core"), col.order],
            paste(args$out_dir, "proactive-core.phenos.tsv", sep="/"),
            col.names=T, row.names=F, sep="\t", quote=F)
write.table(df[which(df$Cohort != "proactive-core"), col.order],
            paste(args$out_dir, "proactive-other.phenos.tsv", sep="/"),
            col.names=T, row.names=F, sep="\t", quote=F)
