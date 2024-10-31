#!/usr/bin/env R

###################################
#            DFCI G2C:            #
# The Germline Genomics of Cancer #
###################################

# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Helper functions for manipulating various kinds of metadata


#' Explode variable by cancer
#'
#' Expands values of a single variable such that there is a one-to-one mapping
#' between cancer diagnoses and the target variable of interest. In other words:
#' this function repeats values of a target variable in situ for individuals with
#' multiple cancer diagnoses
#'
#' @param cancer Character vector with one element per participant, where each
#' element is a semicolon-delimited list of cancer diagnoses for that participant
#' @param x (Optional) vector of values to be exploded
#'
#' @details If no `x` is provided, only the cancer labels will be exploded
#'
#' @returns List of two vectors:
#' * `cancer` : exploded version of input cancer labels, `cancer`
#' * `x` : exploded version of input target vector, `x`
#'
#' @export explode.by.cancer
#' @export
explode.by.cancer <- function(cancer, x=NULL){
  cancer.e <- unlist(strsplit(cancer, split=";", fixed=T))
  x.e <- if(is.null(x)){NULL}else{
    unlist(sapply(1:length(cancer), function(i){
      rep(x[i], length(unlist(strsplit(cancer[i], split=";", fixed=T))))
    }))
  }
  list("cancer" = cancer.e, "x" = x.e)
}


