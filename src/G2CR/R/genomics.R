#!/usr/bin/env R

###################################
#            DFCI G2C:            #
# The Germline Genomics of Cancer #
###################################

# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Various genomics functions


#' Zygosity classification
#'
#' Parse a VCF-style GT character vector and classify the zygosity of each GT
#'
#' @param gt Character vector of GT values formatted according to VCF spec
#'
#' @returns Character vector of zygosity classifications for each value in
#' input `gt`, either `het` or `hom`
#'
#' @details Ignores null/missing alleles and reports `het` if there are two
#' or more unique alleles observed, otherwise reports `hom`
#'
#' @export classify.gt.zygosity
#' @export
classify.gt.zygosity <- function(gt){
  sapply(as.character(gt), function(g){
    n.alleles <- length(unique(setdiff(unlist(strsplit(g, "[/|]")), ".")))
    remap(as.character(n.alleles > 1),
          c("TRUE" = "het", "FALSE" = "hom"))
  })
}
