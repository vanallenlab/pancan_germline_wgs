#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot overall QC summary figures


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
load.constants("all")


##################
# Data Functions #
##################
# Invert dynamic range entries
invert.dynamic.ranges <- function(ss){
  dr.idxs <- grep("dynamic_range", ss$measure)
  ss$value[dr.idxs] <- 1 / ss$value[dr.idxs]
  ss$measure[dr.idxs] <- gsub("dynamic_range", "inverted_dynamic_range", ss$measure[dr.idxs])
  return(ss)
}

# Load summary statistics & organize in sub-dataframes for plotting
load.ss <- function(tsv.in, ref.prefix=NULL, gb.prefixes=c()){
  ss <- read.table(tsv.in, header=F, sep="\t")
  colnames(ss) <- c("analysis", "measure", "value", "n")
  ss <- invert.dynamic.ranges(ss)

  vcs <- c("all", names(var.class.abbrevs))
  ss.l <- lapply(vcs, function(vc){

    # Counts
    if(vc == "all"){
      count.suffix <- "all_variants"
    }else{
      count.suffix <- vc
    }
    prop.rare <- ss[which(ss$analysis == paste("pct_rare", vc, sep=".")), ]
    rare.k <- if(nrow(prop.rare) > 0){round(prop.rare$value * prop.rare$n, 0)}else{0}
    rare.n <- if(nrow(prop.rare) > 0){prop.rare$n}else{0}
    prop.stn <- ss[which(ss$analysis == paste("pct_singletons", vc, sep=".")), ]
    stn.k <- if(nrow(prop.stn) > 0){round(prop.stn$value * prop.stn$n, 0)}else{0}
    stn.n <- if(nrow(prop.stn) > 0){prop.stn$n}else{0}
    count.df <- as.data.frame(rbind(
      ss[which(ss$analysis == paste("site_count", count.suffix, sep=".")), ],
      c(paste("site_count.rare", count.suffix, sep="."), "count", rare.k, rare.n),
      c(paste("site_count.singletons", count.suffix, sep="."), "count", stn.k, stn.n),
      ss[which(ss$analysis == paste("variants_per_genome", vc, sep=".")
               & ss$measure == "median"), ]
    ))

    # Site benchmarking
    sb.df <- as.data.frame(rbind(
      ss[which(ss$analysis == paste(vc, "common_hwe", sep=".")), ],
      ss[which(ss$analysis == paste(vc, "common_ld.any", sep=".")
               & ss$measure == "tag_rate"), ],
      if(vc == "all"){
        data.frame("analysis" = "site_ratios.all", "measure" = NA, "value" = NA, "n" = NA)
      }else{
        ss[which(ss$analysis == paste("site_ratios", vc, sep=".")), ]
      },
      ss[which(ss$analysis == paste(ref.prefix, vc, "common_af_cor", sep=".")), ],
      ss[which(ss$analysis == gsub("^\\.|\\.$", "",
                                   paste(ref.prefix, "sensitivity.common",
                                         if(vc!="all"){vc}, sep="."))), ],
      ss[which(ss$analysis == gsub("^\\.|\\.$", "",
                                   paste(ref.prefix, "ppv.common",
                                         if(vc!="all"){vc}, sep="."))), ]
    ))

    # Genotype benchmarking
    gb.df <- as.data.frame(rbind(
      ss[which(ss$analysis == paste("trio_concordance", vc, sep=".")
               & ss$measure == "median_child_proportion_inherited"), ],
      ss[which(ss$analysis == paste("trio_concordance", vc, sep=".")
               & ss$measure == "child_proportion_inherited_inverted_dynamic_range"), ],
      ss[which(ss$analysis == paste("gt_benchmarking.twin_replicate", vc, sep=".")
               & ss$measure == "median_match_rate"), ],
      ss[which(ss$analysis == paste("gt_benchmarking.twin_replicate", vc, sep=".")
               & ss$measure == "match_rate_inverted_dynamic_range"), ],
      ss[which(ss$analysis == paste("heterozygosity", vc, sep=".")
               & ss$measure == "median"), ],
      ss[which(ss$analysis == paste("heterozygosity", vc, sep=".")
               & ss$measure == "inverted_dynamic_range"), ],
      ss[which(ss$analysis == paste("variants_per_genome", vc, sep=".")
               & ss$measure == "inverted_dynamic_range"), ],
      do.call("rbind", lapply(gb.prefixes, function(gbp){
        rbind(ss[which(ss$analysis %in% paste("gt_benchmarking", gbp, vc, sep=".")
                       & ss$measure == "median_sensitivity"), ],
              ss[which(ss$analysis %in% paste("gt_benchmarking", gbp, vc, sep=".")
                       & ss$measure == "sensitivity_inverted_dynamic_range"), ],
              ss[which(ss$analysis %in% paste("gt_benchmarking", gbp, vc, sep=".")
                       & ss$measure == "median_ppv"), ],
              ss[which(ss$analysis %in% paste("gt_benchmarking", gbp, vc, sep=".")
                       & ss$measure == "ppv_inverted_dynamic_range"), ])
      }))
    ))

    # Inter-class benchmarking
    inter.df <- as.data.frame(rbind(
      do.call("rbind", lapply(c("heterozygosity_cor", "variant_count_cor"), function(suf){
        do.call("rbind", lapply(setdiff(vcs, "all"), function(vc2){
          ordered.vcs <- if(vc != vc2){intersect(vcs, c(vc, vc2))}else{rep(vc, 2)}
          inter.prefix <- paste(ordered.vcs[1], "vs", ordered.vcs[2], sep="_")
          vc.n <- ss$n[which(ss$analysis == paste("variants_per_genome", vc, sep="."))]
          if(length(vc.n) == 0){
            vc.n <- NA
          }
          if(vc == vc2){
            data.frame("analysis" = paste(inter.prefix, suf, sep="."),
                       "measure" = "r2", "value" = 1, "n" = vc.n)
          }else{
            ss[which(ss$analysis == paste(inter.prefix, suf, sep=".")), ]
          }
        }))
      }))
    ))

    # Structured output
    list("counts" = count.df[which(count.df$n > 0 & !is.na(count.df$n)
                                   & !is.infinite(count.df$value)
                                   & !is.na(count.df$value)), ],
         "site_bench" = sb.df[which((sb.df$n > 0 & !is.na(sb.df$n)
                                     & !is.infinite(sb.df$value)
                                     & !is.na(sb.df$value))
                                    | (vc == "all" & sb.df$analysis == "site_ratios.all")), ],
         "gt_bench" = gb.df[which(gb.df$n > 0 & !is.na(gb.df$n)
                                  & !is.infinite(gb.df$value)
                                  & !is.na(gb.df$value)), ],
         "interclass" = inter.df[which(inter.df$n > 0 & !is.na(inter.df$n)
                                       & !is.infinite(inter.df$value)
                                       & !is.na(inter.df$value)), ])
  })
  names(ss.l) <- vcs
  return(ss.l)
}


######################
# Plotting Functions #
######################


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot overall QC summary")
parser$add_argument("--stats", metavar=".tsv", type="character",
                    help="File of all QC summary statistics")
parser$add_argument("--site-ref-prefix", metavar="string", type="character",
                    help="String prefix for site-level benchmarking metrics")
parser$add_argument("--site-ref-title", metavar="string", type="character", default="Ref. cohort",
                    help="Reference dataset title for site-level labels in plots")
parser$add_argument("--sample-benchmarking-prefix", metavar="string",
                    type="character", action="append",
                    help="String prefix for each sample-level benchmarking dataset")
parser$add_argument("--sample-benchmarking-title", metavar="string",
                    type="character", action="append",
                    help="Plot title for each sample-level benchmarking dataset")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV (SINGLE CLASS)
# args <- list("stats" = "~/Downloads/dfci-ufc.sv.v1.InitialVcfQcMetrics.stats/dfci-ufc.sv.v1.InitialVcfQcMetrics.all_qc_summary_metrics.tsv",
#              "site_ref_prefix" = "gnomad-sv_v4.1",
#              "site_ref_title" = "gnomAD-SV v4.1",
#              "sample_benchmarking_prefix" = c("external_srwgs", "external_lrwgs"),
#              "sample_benchmarking_title" = c("External srWGS", "External lrWGS"),
#              "out_prefix" = "~/scratch/ufc_sv.raw.test")

# Load and organize summary stats
ss <- load.ss(args$stats, args$site_ref_prefix, args$sample_benchmarking_prefix)
# TODO: implement this
