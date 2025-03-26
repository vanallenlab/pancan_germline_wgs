#!/usr/bin/env R

###################################
#            DFCI G2C:            #
# The Germline Genomics of Cancer #
###################################

# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Helper functions for customized reading/writing of various files


#' Read & clean .fam file
#'
#' Reads an input .fam file and enforces PLINK standard six-column format
#'
#' @param in.fam Path to input .fam file
#' @param header Does `in.fam` have a header line? \[default: FALSE\]
#'
#' @seealso \url{https://www.cog-genomics.org/plink/1.9/formats#fam}
#'
#' @returns six-column `data.frame`
#'
#' @export load.famfile
#' @export
load.famfile <- function(in.fam, header=FALSE){
  fam <- read.table(in.fam, header=header, sep="\t")
  ncol.input <- ncol(fam)
  fam <- fam[, 1:min(6, ncol.input)]
  fam.colnames <- c("family", "proband", "father", "mother", "sex", "pheno")
  colnames(fam) <- fam.colnames[1:min(6, ncol.input)]
  if(ncol.input < length(fam.colnames)){
    for(k in (ncol.input+1):length(fam.colnames)){
      fam[, fam.colnames[k]] <- 0
    }
  }
  return(fam)
}


