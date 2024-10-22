#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell subroutine to install common dependencies not preinstalled on AoU nodes

# Note that this code is designed to be run inside the AoU Researcher Workbench

# Any positional arguments given to this script will be used to specify which
# subsets of libraries (i.e., which languages) should be installed
if [ $# -eq 0 ]; then
  echo -e "Error: at least one positional argument must be supplied to install_packages.sh"
  echo -e "Current options include: R"
  exit 1
fi


# Install all packages as optioned
for lang in "$@"; do

  case $lang in

    # Install R libraries
    R)
      for lib in argparse beeswarm bedr caret EQL vioplot DescTools; do
        Rscript -e "if(require('$lib') == FALSE){install.packages('$lib', repos='https://cloud.r-project.org')}"
      done

      export rlctools_version=0.1
      Rscript -e "if(require('RLCtools') == TRUE){remove.packages('RLCtools')}; install.packages('~/code/src/RLCtools_$rlctools_version.tar.gz', repos=NULL, type='source')"

      export g2c_version=0.1.0
      Rscript -e "if(require('G2C') == TRUE){remove.packages('G2C')}; install.packages('~/code/src/G2C_$g2c_version.tar.gz', repos=NULL, type='source')"
      ;;

  esac
done