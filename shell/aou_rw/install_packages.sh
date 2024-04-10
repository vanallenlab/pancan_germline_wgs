#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell subroutine to install common dependencies not preinstalled on AoU nodes

# Note that this code is designed to be run inside the AoU Researcher Workbench

# Install R libraries
Rscript -e "install.packages('optparse', repos='https://cloud.r-project.org')"
