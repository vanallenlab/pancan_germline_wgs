#!/usr/bin/env R

###################################
#            DFCI G2C:            #
# The Germline Genomics of Cancer #
###################################

# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Helper functions for plotting


#' Determine point parameters for a scatterplot
#'
#' Automate selection of `pch`, `alpha`, and `cex` for points in a scatterplot
#' depending on the total number of points to be plotted
#'
#' @param n Total number of points to be plotted
#' @param cex.floor Minimum sample size for dynamic `cex` scaling
#' @param cex.ceiling Maximum sample size for dynamic `cex` scaling
#' @param cex.start Largest value of `cex`; used when `n <= cex.floor`
#' @param cex.end Smallest value of `cex`; used when `n >= cex.ceiling`
#' @param pch.cutoff Maximum value of `n` before switching from `pch.big` to `pch.small`
#' @param pch.big Value of `pch` to use when `n <= pch.cutoff`
#' @param pch.small Value of `pch` to use when `n > pch.cutoff`
#' @param alpha.floor Minimum sample size for dynamic `alpha` scaling
#' @param alpha.ceiling Maximum sample size for dynamic `alpha` scaling
#' @param alpha.start Largest value of `alpha`; used when `n <= alpha.floor`
#' @param alpha.end Smallest value of `alpha`; used when `n >= alpha.ceiling`
#'
#' @returns three-element list of `cex`, `pch`, `alpha`
#'
#' @export scatterplot.point.params
#' @export
scatterplot.point.params <- function(n,
                                     cex.floor=1000, cex.ceiling=100000, cex.start=0.6, cex.end=0.05,
                                     pch.cutoff=10000, pch.big=1, pch.small=19,
                                     alpha.floor=100000, alpha.ceiling=1000000, alpha.start=1, alpha.end=0.01){
  # Logscale all counts and ranges
  n <- log10(n)
  cex.range <- sort(log10(c(cex.floor, cex.ceiling)))
  pch.cutoff <- log10(pch.cutoff)
  alpha.range <- sort(log10(c(alpha.floor, alpha.ceiling)))

  # Determine cex
  if(n < cex.range[1]){
    cex <- cex.start
  }else if(n > cex.range[2]){
    cex <- cex.end
  }else{
    cex.pal <- seq(cex.start, cex.end, length.out=101)
    cex.idx <- head(which.min(abs(n - seq(cex.range[1], cex.range[2], length.out=101))), 1)
    cex <- cex.pal[cex.idx]
  }

  # Determine pch
  pch <- if(n <= pch.cutoff){pch.big}else{pch.small}

  # Determine alpha
  if(n < alpha.range[1]){
    alpha <- alpha.start
  }else if(n > alpha.range[2]){
    alpha <- alpha.end
  }else{
    alpha.pal <- seq(alpha.start, alpha.end, length.out=101)
    alpha.idx <- head(which.min(abs(n - seq(alpha.range[1], alpha.range[2], length.out=101))), 1)
    alpha <- alpha.pal[alpha.idx]
  }

  # Return parameters
  return(list("cex" = cex, "pch" = pch, "alpha" = alpha))
}
