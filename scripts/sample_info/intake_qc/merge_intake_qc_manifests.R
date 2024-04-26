#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Merge and harmonize AoU and non-AoU intake QC manifests
# Also assigns G2C IDs to all samples
# Also assigns provisional ancestry labels based on 1kG


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(caret, quietly=TRUE)


##################
# Data Functions #
##################
# Read & clean input .tsv
load.intake.tsv <- function(tsv.in){
  # Read data.frame
  df <- read.table(tsv.in, sep="\t", check.names=F, comment.char="", header=T)

  # Drop unnecessary columns
  cidx.to.drop <- union(grep("^[manta|melt|wham]", colnames(df)),
                        which(colnames(df) %in% c("chrX_count", "chrY_count")))
  if(length(cidx.to.drop) > 0){
    df <- df[, -cidx.to.drop]
  }

  # Remove samples with excessive metadata missingness
  frac.nas <- apply(df, 1, function(vals){length(which(is.na(vals))) / length(vals)})
  ridx.to.drop <- which(frac.nas > 0.5)
  if(length(ridx.to.drop) > 0){
    df <- df[-ridx.to.drop, ]
  }

  return(df)
}


# Merge QC data.frames from AoU and non-AoU samples
merge.qc.dfs <- function(aou.df, other.df, id.prefix="G2C", suffix.length=6){
  # Merge, sort, and deduplicate data
  df <- rbind(aou.df, other.df)
  df <- df[with(df, order(Cohort, Sample)), ]
  df <- df[!duplicated(df), ]

  # Add G2C IDs
  df$G2C_id <- paste(id.prefix, formatC(1:nrow(df), width=suffix.length, flag="0"), sep="")

  return(df)
}


# Infer sex from X/Y ploidy estimates
infer.sex <- function(df, y.tolerance=0.15, x.tolerance=0.25){
  # Begin by estimating copy numbers for X/Y with simple rounding
  n.x <- round(df$chrX_CopyNumber)
  n.y <- round(df$chrY_CopyNumber)

  # Clean up samples with non-binary sex complements
  check.idx <- which(n.x+n.y != 2)
  n.y[check.idx][which(abs(df$chrY_CopyNumber[check.idx] - 1) < 1-y.tolerance)] <- 1
  check.idx.y0 <- intersect(check.idx, which(n.y == 0))
  n.x[check.idx.y0][which(abs(df$chrX_CopyNumber[check.idx.y0] - 2) < 1-x.tolerance)] <- 2
  check.idx.y1 <- intersect(check.idx, which(n.y == 1))
  n.x[check.idx.y1][which(abs(df$chrX_CopyNumber[check.idx.y1] - 1) < 1-x.tolerance)] <- 1

  # Assign sex info to data.frame
  df$inferred_sex <- sapply(paste(n.x, n.y, sep="_"), function(p){
    if(p == "1_1"){"male"}else if(p == "2_0"){"female"}else{"other"}
  })
  df$sex_karyotype <- paste(sapply(n.x, function(n){paste(rep("X", n), collapse="")}),
                            sapply(n.y, function(n){paste(rep("Y", n), collapse="")}),
                            sep="")
  df[which(df$sex_karyotype == "X"), "sex_karyotype"] <- "XO"

  return(df)
}


# Assign provisional ancestries based on KNN clustering vs. HGSV samples
assign.ancestry <- function(df, hgsv.labels.tsv){
  # Read HGSV ground-truth labels
  hgsv.pops <- read.table(hgsv.labels.tsv, header=T, comment.char="",
                          check.names=F, sep="\t")
  hgsv.ids <- as.character(hgsv.pops[, 1])
  hgsv.pops <- as.character(hgsv.pops[, 2])
  names(hgsv.pops) <- hgsv.ids

  # Train KNN on HGSV labels
  train.idx <- which(df$Sample %in% hgsv.ids)
  df.x <- scale(df[, c(grep("^grafpop_GD[1-9]", colnames(df)),
                       which(colnames(df) %in% paste("pct", c("EUR", "AFR", "ASN"), sep="_")))])
  train.x <- df.x[train.idx, ]
  train.y <- hgsv.pops[df$Sample[train.idx]]
  train <- data.frame(train.x, "pop"=train.y)
  set.seed(2024)
  fit <- train(pop ~ ., data=train, method="knn", tuneLength=10)

  # Apply trained KNN to assign ancestry to non-HGSV samples
  pred.pops <- predict(fit, df.x)
  pred.pops[train.idx] <- train.y

  # Add labels to df
  df$intake_qc_pop <- substr(pred.pops, 1, 3)
  df$intake_qc_subpop <- substr(pred.pops, 5, 7)

  return(df)
}


# Clean output data.frame
clean.output <- function(df){
  # Rename columns
  cols.to.rename <- list("Sample" = "original_id",
                         "Cohort" = "cohort",
                         "grafpop_SNPs" = "n_grafpop_snps",
                         "rd_median" = "median_coverage",
                         "rd_mean" = "mean_coverage")
  for(old.name in names(cols.to.rename)){
    colnames(df)[which(colnames(df) == old.name)] <- cols.to.rename[[old.name]]
  }
  cols.to.lower <- c("HQ_HOM", "HQ_HOM_RATE", "HQ_HET", "HQ_HET_RATE",
                     "CHARR", "MEAN_REF_AB_HOM_ALT", "HETEROZYGOSITY_RATE",
                     "INCONSISTENT_AB_HET_RATE")
  to.lower.idx <- which(colnames(df) %in% cols.to.lower)
  colnames(df)[to.lower.idx] <- tolower(colnames(df)[to.lower.idx])
  colnames(df) <- gsub("_CopyNumber$", "_ploidy", colnames(df))

  # Return columns in a specific order
  cols.first <- c("G2C_id", "original_id", "cohort", "inferred_sex", "intake_qc_pop")
  pop.cols <- c(colnames(df)[grep("grafpop", colnames(df))],
                colnames(df)[grep("pct_[A-Z]", colnames(df))],
                "intake_qc_subpop")
  ploidy.cols <- c(colnames(df)[grep("_ploidy", colnames(df))], "sex_karyotype")
  other.cols <- setdiff(colnames(df), c(cols.first, pop.cols, ploidy.cols))
  df[, c(cols.first, pop.cols, ploidy.cols, other.cols)]
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Merge intake QC manifests")
parser$add_argument("--aou-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Intake QC .tsv for AoU samples")
parser$add_argument("--non-aou-tsv", metavar=".tsv", type="character",
                    required=TRUE, help="Intake QC .tsv for non-AoU samples")
parser$add_argument("--hgsv-pop-assignments", metavar=".tsv", type="character",
                    help=paste(".tsv mapping HGSV sample IDs to their ground-",
                               "truth populations. If supplied, all samples ",
                               "will be assigned a provisional ancestry label ",
                               "based on KNN clustering applied to grafpop ",
                               "distances", sep=""))
parser$add_argument("--id-prefix", metavar="string", type="character",
                    default="G2C", help="String prefix for G2C IDs")
parser$add_argument(c("-n", "--suffix-length"), metavar="integer", type="numeric",
                    default=6, help="Length of integer suffixes for G2C IDs")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Path to output .tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("aou_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.tsv.gz",
#              "non_aou_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.tsv.gz",
#              "hgsv_pop_assignments" = "~/scratch/HGSV.ByrskaBishop.sample_populations.tsv",
#              "id_prefix" = "G2C",
#              "suffix_length" = 6,
#              "out_tsv" = "~/scratch/dfci-g2c.intake_qc.merged.test.tsv")

# Load AoU and non-AoU .tsvs
aou.df <- load.intake.tsv(args$aou_tsv)
other.df <- load.intake.tsv(args$non_aou_tsv)

# Merge .tsvs, sort by cohort/ID, and assign G2C IDs
qc.df <- merge.qc.dfs(aou.df, other.df)

# Infer sex from allosome ploidy
qc.df <- infer.sex(qc.df)

# Assign provisional ancestry based on 1kG + grafpop distances, if optioned
if(!is.null(args$hgsv_pop_assignments)){
  qc.df <- assign.ancestry(qc.df, args$hgsv_pop_assignments)
}

# Clean merged data.frame & write to --out-tsv
df.out <- clean.output(qc.df)
colnames(df.out)[1] <- paste("#", colnames(df.out)[1], sep="")
write.table(df.out, args$out_tsv, col.names=T, row.names=F, quote=F, sep="\t")

