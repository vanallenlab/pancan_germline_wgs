#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Curate sample-level phenotype data downloaded from dbGaP for TOPMed BioMe


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2C, quietly=TRUE)
G2C::load.constants("all")


##################
# Data Functions #
##################
# Load & reformat BioMe phenotype data



###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Curate NCI GDC phenotype data")
parser$add_argument("--phenotypes-tsv", metavar=".tsv", type="character",
                    required=TRUE, help="subject_phenotypes.tsv downloaded from dbGaP")
parser$add_argument("--out-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="path to output.tsv")
args <- parser$parse_args()

# # DEV:
# args <- list("phenotypes_tsv" = "~/Desktop/Collins/VanAllen/pancancer_wgs/data_and_cohorts/topmed/biome/phs001644.v2.pht009946.v2.p2.c1.TOPMed_CCDG_BioME_Subject_Phenotypes.HMB-NPU.txt",
#              "out_tsv" = "~/scratch/biome.pheno.dev.tsv")

# Load and clean clinical data
pheno.df <- load.phenotypes(args$phenotypes_tsv)
