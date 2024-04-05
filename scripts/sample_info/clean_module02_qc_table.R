#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Correct WGD info in outputs from GATK-SV module 02
# This is necessary because the WGD scoring script does not have check.names=T


# Set parameters & load libraries
options(scipen=1000, stringsAsFactors=FALSE)

# Only positional argument is QC table from GATK-SV module 02
# This file will be read into memory and updated in situ
tsv.fname <- commandArgs(trailingOnly=TRUE)[1]
df <- read.table(tsv.fname,  header=T, sep="\t", comment.char="", check.names=F)
df$mean_insert_size <- NULL

# Drop VCF outlier columns
# These will be computed later on a cohort-wide basis after accounting for ancestry
df[, grep("outlier", colnames(df))] <- NULL
drop.rows <- which(apply(df[, -c(1)], 1, function(vals){all(vals == "NaN" | is.na(vals))}))
if(length(drop.rows) > 0){
  df <- df[-drop.rows, ]
}

# Find samples with missing WGD scores
bad.wgd.ids <- df[which(df$wgd_score == "NaN" | is.na(df$wgd_score)), 1]

# Correct each sample in serial
for(correct.sid in bad.wgd.ids){
  idx <- which(df[, 1] == correct.sid)
  bad.sid <- gsub("-", ".", correct.sid, fixed=T)
  if(!bad.sid %in% df[, 1]){
    stop(paste("Unable to find duplicate ID (guess: '", bad.sid, 
               "') for sample ", correct.sid, sep=""))
  }else{
    bad.idx <- which(df[, 1] == bad.sid)
    df$wgd_score[idx] <- df$wgd_score[bad.idx]
    df <- df[-c(bad.idx), ]
  }
}

# Update QC .tsv in situ
write.table(df, tsv.fname, col.names=T, row.names=F, quote=F, sep="\t")
