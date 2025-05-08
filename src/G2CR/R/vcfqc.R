#!/usr/bin/env R

###################################
#            DFCI G2C:            #
# The Germline Genomics of Cancer #
###################################

# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Helper functions for VCF QC


#' Clean size or AF break labels
#'
#' Converts size or AF breaks encoded from characters to numeric values
#'
#' @param vals A character vector of break values to be parsed
#'
#' @returns A named vector of numeric values corresponding to `vals`, where
#' element names correspond to the original entries in `vals`
#'
#' @export clean.breaks
#' @export
clean.breaks <- function(vals){
  orig <- vals
  vals <- gsub("[eE]-", " * 10 ^ -", vals)
  vals <- gsub("[eE]\\+", " * 10 ^ ", vals)
  vals <- gsub("[A-Z]|[a-z]|_", "", vals)
  vals <- sapply(vals, function(v){as.numeric(eval(parse(text=v)))})
  names(vals) <- orig
  return(vals)
}


#' Import a compressed summary distribution
#'
#' Import a precomputed compressed summary distribution from a .tsv
#'
#' @param tsv.in Path to input .tsv with compressed distribution
#' @param key.cols Column indexes to be treated as categorical key values \[default: 1-2\]
#'
#' @returns A two-element list of `df`, a dataframe of the distribution encoded
#' in `tsv.in`, and `breaks`, a numeric vector corresponding to the column-axis
#' break values
#'
#' @seealso [G2CR::clean.breaks]
#'
#' @export read.compressed.distrib
#' @export
read.compressed.distrib <- function(tsv.in, key.cols=1:2){
  if(is.null(tsv.in)){return(NULL)}
  df <- read.table(tsv.in, header=T, sep="\t", comment.char="",
                   quote="", check.names=F)
  colnames(df)[1] <- gsub("#", "", colnames(df)[1])
  df[, -key.cols] <- apply(df[, -key.cols], 2, as.numeric)
  dist.breaks <- clean.breaks(colnames(df)[-key.cols])
  list("df" = df, "breaks" = dist.breaks)
}

