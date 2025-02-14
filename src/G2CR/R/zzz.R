#!/usr/bin/env R

###################################
#            DFCI G2C:            #
# The Germline Genomics of Cancer #
###################################

# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Ensure non-standard dependencies are available when package is loaded


.onLoad <- function(libname, pkgname){

  # Set useful global constants
  options(scipen=1000, stringsAsFactors=F, family="sans")

  # Check to make sure RLCtools is available
  if(!require(RLCtools)){
    stop(paste("Dependency `RLCtools` is required for `G2C` but is not found.\n",
              "For more info, see https://github.com/RCollins13/RLCtools\n", sep=""))
  }
}
