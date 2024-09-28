#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data for cohorts downloaded from NIH/NCI GDC


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2C, quietly=TRUE)
G2C::load.constants("all")

# Declare constants used in variable parsing
race.map <- c("american indian or alaska native" = "native_american",
              "black or african american" = "black",
              "not reported" = NA,
              "Unknown" = NA)
eth.map <- c("hispanic or latino" = "hispanic",
             "not hispanic or latino" = "not_hispanic",
             "not reported" = NA,
             "Unknown" = NA)
vital.map <- c("Alive" = 1,
               "Dead" = 0,
               "Unknown" = NA)
origin.map <- c("prostate_gland" = "prostate",
                "corpus_uteri" = "uterus",
                "endometrium" = "uterus",
                "larynx" = "oral_cavity",
                "brain" = "cns",
                "unknown" = "other",
                "lip" = "oral_cavity",
                "oropharynx" = "oral_cavity",
                "gum" = "oral_cavity",
                "cheek_mucosa" = "oral_cavity",
                "head_face_or_neck" = "oral_cavity",
                "tonsil" = "oral_cavity",
                "skin" = "melanoma",
                "colon" = "colorectal",
                "rectum" = "colorectal",
                "ampulla_of_vater" = "other",
                "small_intestine" = "other",
                "gallbladder" = "other",
                "not_reported" = "other")
origin.grep.map <- c("tongue" = "oral_cavity",
                     "mouth" = "oral_cavity",
                     "bone" = "other",
                     "connective" = "other",
                     "soft_tissue" = "other",
                     "bile_duct" = "other",
                     "rectosigmoid" = "colorectal",
                     "nasal" = "oral_cavity",
                     "oral_cavity" = "oral_cavity")
dx.override.map <- c("glioblastoma" = "cns")
stage.map <- c("Unknown" = "unknown",
               "Not Reported" = "unknown",
               "NA" = "unknown",
               "Not R" = "unknown",
               "X" = "unknown",
               "0" = "0",
               "1" = "I",
               "2" = "II",
               "3" = "III",
               "4" = "IV")
met.map <- c("Distant Metastasis" = TRUE,
             "Metastasis, NOS" = TRUE,
             "No Metastasis" = FALSE,
             "Unknown" = FALSE)
grade.map <- c("Not Reported" = "unknown",
               "X" = "unknown",
               "B" = "1",
               "High Grade" = "4",
               "Unknown" = "unknown",
               "NA" = "unknown",
               "3;GX" = "3")
specimen.map <- c("Blood Derived Normal" = "blood",
                  "Buccal Cell Normal" = "buccal_cells",
                  "Solid Tissue Normal" = "solid_tissue")


##################
# Data Functions #
##################
# Simplify tissue or organ of origin
simplify.origin <- function(vals){
  vals <- remap(vals, origin.map)
  sapply(vals, function(o){
    for(organ in c("lung", "pancreas", names(origin.grep.map))){
      if(grepl(paste("(^|_)", organ, "(_|$)", sep=""), o)){
        return(remap(organ, origin.grep.map))
      }
    }
    return(o)
  })
}

# "Squash" rows in a data frame
squash.duplicates <- function(dup.df){
  # Start by handling all columns in a generic manner
  dedup.vals <- apply(dup.df, 2, function(vals){
    vals <- unlist(vals)
    if(all(is.na(vals))){
      NA
    }else{
      vals <- unique(vals[which(!is.na(vals))])
      if(length(vals) == 1){
        vals[1]
      }else{
        if(all(!is.na(as.numeric(vals)))){
          median(as.numeric(vals))
        }else{
          paste(sort(vals), collapse=";")
        }
      }
    }
  })

  # Update select columns with more nuanced logic
  min.cols <- c("age_at_diagnosis", "year_of_birth", "days_to_birth")
  max.cols <- c("days_to_death", "days_to_last_follow_up")
  for(col in c(min.cols, max.cols)){
    cv <- as.numeric(unlist(dup.df[, col]))
    if(all(is.na(cv))){
      NA
    }else{
      if(col %in% min.cols){
        dedup.vals[col] <- min(cv, na.rm=T)
      }else{
        dedup.vals[col] <- max(cv, na.rm=T)
      }
    }
  }
  vital <- as.character(unlist(dup.df$vital_status))
  if(any(!is.na(vital))){
    vital <- vital[which(!is.na(vital))]
    dedup.vals$vital_status <- if("Dead" %in% vital){"Dead"}else{"Alive"}
  }

  as.data.frame(t(dedup.vals))
}

# Load & curate clinical.tsv from NCI GDC
load.clinical.data <- function(tsv.in){
  # Read data
  df <- read.table(tsv.in, sep="\t", comment.char="", quote="", header=T,
                   na.strings="'--")

  # Only keep relevant columns
  keep.cols <- c("case_submitter_id", "gender", "race", "ethnicity",
                 "age_at_diagnosis", "year_of_birth", "vital_status",
                 "days_to_birth", "days_to_death", "days_to_last_follow_up",
                 "tissue_or_organ_of_origin", "primary_diagnosis", "icd_10_code",
                 "ajcc_pathologic_stage", "ajcc_clinical_stage",
                 "ajcc_pathologic_t", "gleason_grade_group", "tumor_grade",
                 "site_of_resection_or_biopsy", "primary_disease",
                 "metastasis_at_diagnosis", "metastasis_at_diagnosis_site")
  df <- df[, keep.cols]

  # Deduplicate based on columns above
  df <- df[!duplicated(df), ]

  # Check for duplicate sample IDs, and deduplicate further if needed
  dup.ids <- unique(df$case_submitter_id[which(duplicated(df$case_submitter_id))])
  if(length(dup.ids) > 0){
    for(id in dup.ids){
      # Extract original rows for duplicate ID from main data frame
      dup.ridx <- which(df$case_submitter_id == id)
      dup.df <- df[dup.ridx, ]
      df <- df[-dup.ridx, ]

      # Integrate information across multiple duplicate rows
      dedup.df <- squash.duplicates(dup.df)
      df <- as.data.frame(rbind(df, dedup.df))
    }
  }

  # Enforce column types prior to parsing
  # This is necessary because sometimes the deduplication process
  # causes columns to be converted to lists
  str.cols <- c("case_submitter_id", "gender", "race", "ethnicity",
                "vital_status", "tissue_or_organ_of_origin", "primary_diagnosis",
                "icd_10_code", "ajcc_pathologic_stage", "ajcc_clinical_stage",
                "ajcc_pathologic_t", "gleason_grade_group", "tumor_grade",
                "site_of_resection_or_biopsy", "primary_disease",
                "metastasis_at_diagnosis", "metastasis_at_diagnosis_site")
  df[, str.cols] <- apply(df[, str.cols], 2, function(v){as.character(unlist(v))})
  num.cols <- c("age_at_diagnosis", "year_of_birth", "days_to_birth",
                "days_to_death", "days_to_last_follow_up")
  df[, num.cols] <- apply(df[, num.cols], 2, function(v){as.numeric(unlist(v))})

  # Reassign simple columns requiring little-to-no transformation
  df$Sample <- df$case_submitter_id
  df$reported_sex <- tolower(df$gender)
  na.sex.idx <- which(!df$reported_sex %in% c("male", "female", "other"))
  df$reported_sex[na.sex.idx] <- NA
  df$birth_year <- as.numeric(df$year_of_birth)
  df$vital_status <- as.numeric(remap(df$vital_status, vital.map))
  df$years_to_last_contact <- as.numeric(df$days_to_last_follow_up) / 365
  df$years_to_last_contact[which(df$years_to_last_contact < 0)] <- 0
  df$years_left_censored <- 0
  df$years_left_censored[which(is.na(df$years_to_last_contact))] <- NA
  df$cancer_icd10 <- df$icd_10_code

  # Parse age with variable units
  df$age <- as.numeric(df$age_at_diagnosis)
  age.in.days <- df$age > 150
  if(any(age.in.days)){
    df$age[which(age.in.days)] <- round(df$age[which(age.in.days)] / 365, 1)
  }

  # Compute age at last contact
  df$age_at_last_contact <- round(df$age + (as.numeric(df$days_to_last_follow_up) / 365), 1)

  # Assign race + ethnicity
  df$race <- remap(df$race, race.map)
  df$ethnicity <- remap(df$ethnicity, eth.map)
  if(all(is.na(df$ethnicity))){
    df$rroe <- df$race
  }else{
    df$rroe <- apply(df[, c("race", "ethnicity")], 1, paste, collapse="_")
  }
  df$reported_race_or_ethnicity <-
    gsub("^NA_", "", gsub("_NA$", "", gsub("[\ ]+", "_", df$rroe)))
  rroe.nas <- is.na(df$reported_race_or_ethnicity)
  if(any(rroe.nas)){
    df$reported_race_or_ethnicity[which(rroe.nas)] <- "unknown"
  }
  no.race <- is.na(df$race)
  if(any(no.race)){
    df$reported_race_or_ethnicity[which(no.race)] <- "unknown"
  }

  # Curate cancer diagnosis information
  dx.cols <- c("tissue_or_organ_of_origin", "primary_diagnosis",
               "primary_disease", "ajcc_pathologic_stage")
  df[, setdiff(dx.cols, "ajcc_pathologic_stage")] <-
    apply(df[, setdiff(dx.cols, "ajcc_pathologic_stage")], 2, function(v){
      gsub("_$", "", gsub("[,| ]+", "_", gsub(", nos;", ";", gsub(", nos$", "", tolower(v)))))
    })
  simple_origin <- sapply(df$tissue_or_organ_of_origin, function(v){
    paste(setdiff(unlist(strsplit(v, split=";")), c("unknown", "not_reported")),
          collapse=";")
  })
  df$cancer <- simplify.origin(simple_origin)
  primary.dx.remap <- df$primary_diagnosis %in% names(dx.override.map)
  if(any(primary.dx.remap)){
    p.dx.idx <- which(primary.dx.remap)
    df$cancer[p.dx.idx] <- remap(df$primary_diagnosis[p.dx.idx], dx.override.map)
  }
  missing.dx <- as.logical(as.numeric(apply(df[, dx.cols], 1, function(v){all(is.na(v))}))
                           + as.numeric(df$cancer == ""))
  if(any(missing.dx)){
    df$cancer[which(missing.dx)] <- "other"
  }
  df$odx <- apply(df[, dx.cols], 1, paste, collapse="_")
  df$odx <- gsub("(^|_)NA(_|$)", "_", gsub("[\ ]+", "_", df$odx))
  df$original_dx <- tolower(gsub("(^|_)NA(_|$)", "", df$odx))

  # Curate stage, defaulting to AJCC T category if overall stage is not provided
  na.stage <- is.na(df$ajcc_pathologic_stage)
  if(any(na.stage)){
    df$ajcc_pathologic_stage[which(na.stage)] <- df$ajcc_clinical_stage[which(na.stage)]
  }
  df$stage <- gsub("[ABC][1-9]*$", "", gsub("^Stage ", "", df$ajcc_pathologic_stage))
  na.stage <- is.na(df$stage)
  if(any(na.stage)){
    df$stage[which(na.stage)] <-
      gsub("^T", "", gsub("[a-z]*$", "", df$ajcc_pathologic_t[which(na.stage)]))
  }
  df$stage <- remap(as.character(df$stage), stage.map)

  # Add metastatic indicator
  df$metastatic <- apply(df[, dx.cols], 1, function(v){any(grepl("metastat", v))})
  df$metastatic[which(df$stage == "IV")] <- TRUE
  df$metastasis_at_diagnosis <- remap(df$metastasis_at_diagnosis, met.map)
  df$metastatic[which(as.logical(df$metastasis_at_diagnosis))] <- TRUE
  df$metastatic <- as.integer(df$metastatic)

  # Parse grade
  df$grade <- remap(gsub("^G", "", df$tumor_grade), grade.map)

  # Return formatted columns
  df[, c("Sample", "reported_sex", "reported_race_or_ethnicity", "age",
         "birth_year", "vital_status", "age_at_last_contact",
         "years_to_last_contact", "years_left_censored", "cancer", "stage",
         "metastatic", "grade", "cancer_icd10", "original_dx")]
}

# Assign WGS tissue to each sample, where possible
get.wgs.tissue <- function(df, sample.tsv){
  # Update tissue labels if --biospecimen-tsv is provided
  if(!is.null(sample.tsv)){
    # Load specimen information, subset to normals, and deduplicate
    s.df <- read.table(sample.tsv, header=T, sep="\t", comment.char="", quote="",
                       na.strings="'--")
    keep.cols <- c("case_submitter_id", "preservation_method", "sample_type",
                   "specimen_type", "tissue_type")
    norm.idx <- union(grep("normal", s.df$sample_type, ignore.case=TRUE),
                      grep("normal", s.df$tissue_type, ignore.case=TRUE))
    keep.rows <- intersect(norm.idx, which(s.df$case_submitter_id %in% df$Sample))
    s.df <- s.df[keep.rows, keep.cols]
    s.df <- s.df[c(which(s.df$preservation_method != "Unknown"),
                   which(s.df$preservation_method == "Unknown")), ]
    s.df <- s.df[c(grep("blood", s.df$sample_type, ignore.case=TRUE),
                   grep("blood", s.df$sample_type, ignore.case=TRUE, invert=TRUE)), ]
    s.df <- s.df[which(!duplicated(s.df$case_submitter_id)), ]
    rownames(s.df) <- s.df$case_submitter_id

    # Clean columns & build DNA source map
    s.df$preservation_method <- gsub("[ ]+", "_", tolower(s.df$preservation_method))
    s.df$preservation_method[which(s.df$preservation_method == "unknown")] <- NA
    s.df$sample_type <- remap(s.df$sample_type, specimen.map)
    dna.map <- apply(s.df[, c("preservation_method", "sample_type")], 1, paste, collapse="_")
    dna.map <- gsub("^NA_", "", dna.map[which(!is.na(dna.map))])
  }else{
    dna.map <- c()
  }

  # Map DNA source where possible
  df$wgs_tissue <- remap(df$Sample, dna.map, default.value="unknown")
  df$wgs_tissue[which(is.na(df$wgs_tissue))] <- "unknown"

  return(df)
}

# Load & curate exposure/biometric data
load.exposure.data <- function(tsv.in){
  # Read data
  df <- read.table(tsv.in, sep="\t", comment.char="", quote="", header=T,
                   na.strings="'--")
  if(nrow(df) == 0){
    return(NULL)
  }

  # Only keep relevant columns
  keep.cols <- c("case_submitter_id", "bmi", "height", "weight",
                 "cigarettes_per_day", "pack_years_smoked",
                 "tobacco_smoking_onset_year", "tobacco_smoking_quit_year",
                 "years_smoked", "tobacco_smoking_status")
  df <- df[, keep.cols]

  # Deduplicate based on columns above
  df <- df[!duplicated(df), ]

  # Reassign case ID
  df$Sample <- df$case_submitter_id

  # Infer smoking history
  smoke.cols <- c("cigarettes_per_day", "pack_years_smoked",
                  "tobacco_smoking_onset_year", "tobacco_smoking_quit_year",
                  "years_smoked")
  no.smoke.data <- apply(df[, smoke.cols], 1, function(v){all(is.na(v))})
  df$smoking_history <- as.integer(apply(df[, smoke.cols], 1, function(v){
    sum(as.numeric(v), na.rm=T) > 0
    }))
  if(any(no.smoke.data)){
    df$smoking_history[which(no.smoke.data)] <- NA
  }
  non.smoker <- grepl("lifelong non-smoker", tolower(df$tobacco_smoking_status))
  if(any(non.smoker)){
    df$smoking_history[which(non.smoker)] <- 0
  }
  yes.smoker <- grepl("reformed smoker|current smoker", tolower(df$tobacco_smoking_status))
  if(any(yes.smoker)){
    df$smoking_history[which(yes.smoker)] <- 1
  }

  # Return exposure data
  df[, c("Sample", "height", "weight", "bmi", "smoking_history")]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate NCI GDC phenotype data")
parser$add_argument("--clinical-tsv", metavar=".tsv", type="character",
                    required=TRUE, help="clinical.tsv downloaded from NCI GDC")
parser$add_argument("--exposure-tsv", metavar=".tsv", type="character",
                    help="exposure.tsv downloaded from NCI GDC")
parser$add_argument("--biospecimen-tsv", metavar=".tsv", type="character",
                    help="sample.tsv downloaded from NCI GDC")
parser$add_argument("--cohort", metavar="string", type="character",
                    help=paste("Optional cohort specifier; will be added as a",
                               "column to --out-tsv"))
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("clinical_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/cptac/clinical.cases_selection.2024-03-27/clinical.tsv",
#              "exposure_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/cptac/clinical.cases_selection.2024-03-27/exposure.tsv",
#              "biospecimen_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/cptac/biospecimen.project-cptac-3.2024-08-20/sample.tsv",
#              "cohort" = "cptac",
#              "out_tsv" = "~/scratch/cptac.pheno.dev.tsv")

# Load and clean clinical data
clin.df <- load.clinical.data(args$clinical_tsv)

# Add DNA source, if available
clin.df <- get.wgs.tissue(clin.df, args$biospecimen_tsv)

# Load and clean exposure data
exp.df <- load.exposure.data(args$exposure_tsv)

# Merge dataframes and write to --out-tsv
if(!is.null(exp.df)){
  out.df <- merge(clin.df, exp.df, by="Sample", all.x=T, all.y=F, sort=F)
}else{
  out.df <- clin.df
  out.df[, c("height", "weight", "bmi", "smoking_history")] <- NA
}
if(!is.null(args$cohort)){
  out.df$Cohort <- args$cohort
}
out.cols <- c("Sample", "Cohort", "reported_sex", "reported_race_or_ethnicity",
              "age", "birth_year", "vital_status", "age_at_last_contact",
              "years_to_last_contact", "years_left_censored", "height",
              "weight", "bmi", "cancer", "stage", "metastatic", "grade",
              "smoking_history", "cancer_icd10", "original_dx", "wgs_tissue")
write.table(out.df[order(out.df$Sample), intersect(out.cols, colnames(out.df))],
            args$out_tsv, sep="\t", col.names=T, row.names=F, quote=F)
