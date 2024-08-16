#!/usr/bin/env R

###################################
#            DFCI G2C:            #
# The Germline Genomics of Cancer #
###################################

# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


# Project-wide constants


#' Load Study Constants
#'
#' Load a subset of constants used throughout G2C
#'
#' @param susbet Vector of constant groups to load See `Details` for options.
#' @param envir Environment passed to [base::assign] \[default: .GlobalEnv\]
#'
#' @details Recognized values for `subset` include:
#' * `colors` : all color palettes used
#' * `scales` : all scales and scale labels
#' * `names` : names of various variables
#' * `all` : load all constants
#'
#' @examples
#' # Load list of color palettes
#' get.constants("colors");
#'
#' # Load scales and names colors
#' get.constants(c("scales", "names"))
#'
#' @seealso [base::assign]
#'
#' @export load.constants
#' @export
load.constants <- function(subset, envir=.GlobalEnv){
  # Define colors
  # TODO: add cancer type colors
  male.colors <- c("dark2" = "#2A5869",
                   "dark1" = "#3F839D",
                   "main" = "#54AFD1",
                   "light1" = "#76BFDA",
                   "light2" = "#98CFE3")
  female.colors <- c("dark2" = "#73475E",
                     "dark1" = "#AD6B8C",
                     "main" = "#E68EBB",
                     "light1" = "#EBA5C9",
                     "light2" = "#F0BBD6")
  AFR.colors <- c("dark3" = "#796624",
                  "dark2" = "#A18730",
                  "dark1" = "#C9A93C",
                  "main" = "#F1CB48",
                  "light1" = "#F3D160",
                  "light2" = "#F5D878",
                  "light3" = "#F7E091")
  AMR.colors <- c("dark2" = "#652223",
                  "dark1" = "#973234",
                  "main" = "#C94345",
                  "light1" = "#D4696A",
                  "light2" = "#DF8E8F")
  EAS.colors <- c("dark2" = "#5B6519",
                  "dark1" = "#889826",
                  "main" = "#B5CA32",
                  "light1" = "#C4D55B",
                  "light2" = "#D3DF84")
  EUR.colors <- c("dark2" = "#4E6774",
                  "dark1" = "#749AAE",
                  "main" = "#9BCDE8",
                  "light1" = "#AFD7ED",
                  "light2" = "#C3E1F1")
  SAS.colors <- c("dark2" = "#4D2D4E",
                  "dark1" = "#734474",
                  "main" = "#995A9B",
                  "light1" = "#AD7BAF",
                  "light2" = "#C29CC3")
  DEL.colors <- c("dark2" = "#6A1D13",
                  "dark1" = "#9F2B1C",
                  "main" = "#D43925",
                  "light1" = "#DD6151",
                  "light2" = "#EEB0A8")
  DUP.colors <- c("dark2" = "#123B59",
                  "dark1" = "#1A5985",
                  "main" = "#2376B2",
                  "light1" = "#4F91C1",
                  "light2" = "#A7C8E0")
  CNV.colors <- c("dark2" = "#3A2D59",
                  "dark1" = "#574385",
                  "main" = "#7459B2",
                  "light1" = "#907AC1",
                  "light2" = "#C7BDE0")
  INS.colors <- c("dark2" = "#6A3A70",
                  "dark1" = "#9F57A8",
                  "main" = "#D474E0",
                  "light1" = "#DD90E6",
                  "light2" = "#EEC7F3")
  INV.colors <- c("dark2" = "#7D4A0F",
                  "dark1" = "#BB6E16",
                  "main" = "#FA931E",
                  "light1" = "#FBA94B",
                  "light2" = "#FDD4A5")
  CPX.colors <- c("dark2" = "#397246",
                  "dark1" = "#55AA69",
                  "main" = "#71E38C",
                  "light1" = "#8DE9A3",
                  "light2" = "#C6F4D1")
  CTX.colors <- c("dark2" = "#1D3923",
                  "dark1" = "#2B5635",
                  "main" = "#397246",
                  "light1" = "#618E6B",
                  "light2" = "#B0C7B5")
  BND.colors <- CTX.colors
  csq.colors <- c("synonymous" = "#AAAAAA",
                  "missense" = "#FF6103",
                  "lof" = "#9D1309")
  DFCI.colors <- c("darkblue" = "#003354",
                   "midblue" = "#02679A",
                   "lightblue" = "#3AC6F3",
                   "paleblue" = "#D4EEFF",
                   "yellow" = "#F89820")
  colors <- list(
    "sex.colors" = c("MALE" = male.colors[["main"]],
                     "FEMALE" = female.colors[["main"]],
                     "OTHER" = "#666245"),
    "MALE.colors" = male.colors,
    "FEMALE.colors" = female.colors,
    "pop.colors" = c("AFR" = AFR.colors[["main"]],
                     "AMR" = AMR.colors[["main"]],
                     "EAS" = EAS.colors[["main"]],
                     "EUR" = EUR.colors[["main"]],
                     "SAS" = SAS.colors[["main"]],
                     "OTH" = "gray"),
    "pop.palettes" = list("AFR" = AFR.colors,
                          "AMR" = AMR.colors,
                          "EAS" = EAS.colors,
                          "EUR" = EUR.colors,
                          "SAS" = SAS.colors,
                          "OTH" = "gray"),
    "AFR.colors" = AFR.colors,
    "AMR.colors" = AMR.colors,
    "EAS.colors" = EAS.colors,
    "EUR.colors" = EUR.colors,
    "SAS.colors" = SAS.colors,
    "subpop.colors" = c("ACB" = AFR.colors[["light1"]],
                        "ASW" = AFR.colors[["main"]],
                        "ESN" = AFR.colors[["dark3"]],
                        "GWD" = AFR.colors[["light2"]],
                        "LWK" = AFR.colors[["dark2"]],
                        "MSL" = AFR.colors[["light3"]],
                        "YRI" = AFR.colors[["dark1"]],
                        "CLM" = AMR.colors[["light2"]],
                        "MXL" = AMR.colors[["light1"]],
                        "PEL" = AMR.colors[["dark1"]],
                        "PUR" = AMR.colors[["main"]],
                        "CDX" = EAS.colors[["light2"]],
                        "CHB" = EAS.colors[["dark2"]],
                        "CHS" = EAS.colors[["main"]],
                        "JPT" = EAS.colors[["dark1"]],
                        "KHV" = EAS.colors[["light1"]],
                        "CEU" = EUR.colors[["main"]],
                        "FIN" = EUR.colors[["light2"]],
                        "GBR" = EUR.colors[["dark2"]],
                        "IBS" = EUR.colors[["dark1"]],
                        "TSI" = EUR.colors[["light1"]],
                        "BEB" = SAS.colors[["dark1"]],
                        "GIH" = SAS.colors[["light1"]],
                        "ITU" = SAS.colors[["dark2"]],
                        "PJL" = SAS.colors[["main"]],
                        "STU" = SAS.colors[["light2"]]),
    "DEL.colors" = DEL.colors,
    "DUP.colors" = DUP.colors,
    "CNV.colors" = CNV.colors,
    "INS.colors" = INS.colors,
    "INV.colors" = INV.colors,
    "CPX.colors" = CPX.colors,
    "BND.colors" = BND.colors,
    "CTX.colors" = CTX.colors,
    "OTH.colors" = BND.colors,
    "sv.colors" = c("DEL" = DEL.colors[["main"]],
                    "DUP" = DUP.colors[["main"]],
                    "CNV" = CNV.colors[["main"]],
                    "INS" = INS.colors[["main"]],
                    "INV" = INV.colors[["main"]],
                    "CPX" = CPX.colors[["main"]],
                    "BND" = BND.colors[["main"]],
                    "CTX" = CTX.colors[["main"]],
                    "OTH" = BND.colors[["main"]]),
    "sv.palettes" = list("DEL" = DEL.colors,
                         "DUP" = DUP.colors,
                         "CNV" = CNV.colors,
                         "INS" = INS.colors,
                         "INV" = INV.colors,
                         "CPX" = CPX.colors,
                         "BND" = BND.colors,
                         "CTX" = CTX.colors,
                         "OTH" = BND.colors),
    "csq.colors" = csq.colors,
    "DFCI.colors" = DFCI.colors,
    "stage.colors" = c("0" = "white",
                       "1" = "#F8FAA7",
                       "2" = "#FFCC66",
                       "3" = "#FE8002",
                       "4" = "#ED3823"),
    "relative.colors" = c("duplicates" = "#9D1309",
                          "parent-child" = "#FF6103",
                          "siblings" = "#FFB14D",
                          "unrelated" = "#AAAAAA"))

  # Define scales
  logscale.major <- 10^(-10:10)
  contig.lengths <- c("chr1" = 248956422,
                      "chr2" = 242193529,
                      "chr3" = 198295559,
                      "chr4" = 190214555,
                      "chr5" = 181538259,
                      "chr6" = 170805979,
                      "chr7" = 159345973,
                      "chr8" = 145138636,
                      "chr9" = 138394717,
                      "chr10" = 133797422,
                      "chr11" = 135086622,
                      "chr12" = 133275309,
                      "chr13" = 114364328,
                      "chr14" = 107043718,
                      "chr15" = 101991189,
                      "chr16" = 90338345,
                      "chr17" = 83257441,
                      "chr18" = 80373285,
                      "chr19" = 58617616,
                      "chr20" = 64444167,
                      "chr21" = 46709983,
                      "chr22" = 50818468,
                      "chrX" = 156040895,
                      "chrY" = 57227415)
  scales <- list(
    "logscale.major" = logscale.major,
    "logscale.major.bp" = 10^(0:9),
    "logscale.major.bp.labels" = c(sapply(c("bp", "kb", "Mb"),
                                          function(suf){paste(c(1, 10, 100), suf, sep="")}),
                                   "1 Gb"),
    "logscale.demi" = as.numeric(sapply(logscale.major, function(e){c(1, 5)*e})),
    "logscale.demi.bp" = as.numeric(sapply(10^(0:9), function(e){c(1, 5)*e})),
    "logscale.demi.bp.labels" = c(paste(c(1, 5, 10, 50, 100, 500), "bp", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "kb", sep=""),
                                  paste(c(1, 5, 10, 50, 100, 500), "Mb", sep=""),
                                  paste(c(1, 5), "Gb", sep="")),
    "logscale.minor" = as.numeric(sapply(logscale.major, function(e){(1:9)*e})),
    "yearscale.major" = 0:100 * 365,
    "yearscale.demi" = seq(0, 100, 0.5) * 365,
    "yearscale.minor" = seq(0, 100, 1/12) * 365,
    "contig.lengths" = contig.lengths
  )

  # Define names
  all.names <- list(
    "pop.abbreviations" = c("AFR" = "Afr.",
                            "AMR" = "Amer.",
                            "EAS" = "E. Asn.",
                            "EUR" = "Eur.",
                            "SAS" = "S. Asn.",
                            "OTH" = "Other"),
    "pop.names.short" = c("AFR" = "African",
                          "AMR" = "American",
                          "EAS" = "E. Asian",
                          "EUR" = "European",
                          "SAS" = "S. Asian",
                          "OTH" = "Other"),
    "pop.names.long" = c("AFR" = "African/African-American",
                         "AMR" = "Latino/admixed American",
                         "EAS" = "East Asian",
                         "EUR" = "European",
                         "SAS" = "South Asian",
                         "OTH" = "Other/unknown"),
    "sex.names" = c("MALE" = "Male",
                    "FEMALE" = "Female",
                    "OTHER" = "Other"),
    "stage.names" = c("0" = "",
                      "1" = "I",
                      "2" = "II",
                      "3" = "III",
                      "4" = "IV"),
    "relative.names" = c("duplicates" = "Identical",
                         "parent-child" = "Parent-child",
                         "siblings" = "Siblings",
                         "unrelated" = "Unrelated"),
    "sv.abbreviations" = c("DEL" = "Del.",
                           "DUP" = "Dup.",
                           "CNV" = "mCNV",
                           "INS" = "Ins.",
                           "INV" = "Inv.",
                           "CPX" = "Complex",
                           "CTX" = "Tloc.",
                           "OTH" = "Other"),
    "sv.names" = c("DEL" = "Deletion",
                   "DUP" = "Duplication",
                   "CNV" = "mCNV",
                   "INS" = "Insertion",
                   "INV" = "Inversion",
                   "CPX" = "Complex SV",
                   "CTX" = "Translocation",
                   "OTH" = "Other SV"),
    "csq.names.short" = c("synonymous" = "Syn.",
                          "missense" = "Mis.",
                          "lof" = "LoF")
  )

  # Assign constants to global environment
  if(length(intersect(subset, c("colors", "all"))) > 0){
    for(variable in names(colors)){
      assign(variable, colors[[variable]], envir=envir)
    }
  }
  if(length(intersect(subset, c("scales", "all"))) > 0){
    for(variable in names(scales)){
      assign(variable, scales[[variable]], envir=envir)
    }
  }
  if(length(intersect(subset, c("names", "all"))) > 0){
    for(variable in names(all.names)){
      assign(variable, all.names[[variable]], envir=envir)
    }
  }
}
