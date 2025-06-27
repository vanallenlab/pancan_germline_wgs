#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Benchmark sample-/genotype-level variant overlap between two callsets


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)


##################
# Data functions #
##################
# Read variant map and label with vc, vsc, and allele frequency designation
load.vid.map <- function(vid.map.in, sites.bed.in, common.af=0.01){
  vid.map <- read.table(vid.map.in, header=F, sep="\t")[, 1:2]
  colnames(vid.map) <- c("source_vid", "target_vid")
  bed <- read.table(sites.bed.in, header=T, sep="\t", comment.char="")
  bed$freq <- remap(as.character(bed$af < common.af),
                    c("TRUE" = paste("lt", format(common.af, scientific=T), sep=""),
                      "FALSE" = paste("ge", format(common.af, scientific=T), sep="")))
  merge(vid.map, bed[, c("vid", "class", "subclass", "freq")],
        by.x="source_vid", by.y="vid", all.x=T, all.y=F, sort=F)
}

# Load a map of source & target sample IDs
load.sid.map <- function(sid.map.in){
  sid.map <- read.table(sid.map.in, header=F, sep="\t")
  colnames(sid.map) <- c("source_sid", "target_sid")
  return(sid.map)
}

# Load genotypes for a single sample
load.gts <- function(gt.tsv.in, prefix, elig.vids=NULL){
  gts <- tryCatch(read.table(gt.tsv.in, header=F, sep="\t"),
                  error=function(e){NULL},
                  warning=function(w){NULL})
  if(is.null(gts)){
    return(NULL)
  }
  colnames(gts) <- paste(prefix, c("vid", "gt"), sep="_")
  if(!is.null(elig.vids)){
    gts <- gts[which(gts[, 1] %in% elig.vids), ]
  }
  gts[which(!gts[, 1] %in% c("0", "0/0", ".", "./.", "./0")), ]
}

# Benchmark source & target genotypes
benchmark.gts <- function(source.gts, target.gts, vid.map, report.by.gt=FALSE){
  # Compute concordance for all genotypes
  hits <- merge(merge(vid.map, source.gts, all=F, sort=F),
                target.gts, all.x=T, all.y=F, sort=F)
  hits$zygosity <- if(report.by.gt){hits$source_gt}else{classify.gt.zygosity(hits$source_gt)}
  hits$match <- "none"
  hits$match[which(!is.na(hits$target_gt))] <- "carrier"
  hits$match[which(hits$source_gt == hits$target_gt)] <- "gt"

  # Enumerate GT concordance rate per (vc, vsc, freq)
  do.call("rbind", lapply(sort(unique(hits$class)), function(vc){
    vc.df <- hits[which(hits$class == vc), ]

    vc.res <- do.call("rbind", lapply(sort(unique(vc.df$subclass)), function(vsc){
      vsc.df <- vc.df[which(vc.df$subclass == vsc), ]

      vsc.res <- do.call("rbind", lapply(sort(unique(vsc.df$freq)), function(freq){
        freq.df <- vsc.df[which(vsc.df$freq == freq), ]

        freq.res <- do.call("rbind", lapply(sort(unique(freq.df$zygosity)), function(zyg){
          zyg.df <- freq.df[which(freq.df$zygosity == zyg), ]

          c(zyg, sapply(c("none", "carrier", "gt"), function(mt){
            sum(zyg.df$match == mt)
          }))

        }))
        cbind(rep(freq, nrow(freq.res)), freq.res)

      }))
      cbind(rep(vsc, nrow(vsc.res)), vsc.res)

    }))

    cbind(rep(vc, nrow(vc.res)), vc.res)
  }))
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Benchmark sample genotypes between callsets")
parser$add_argument("--variant-map", metavar=".tsv", required=TRUE,
                    help=paste(".tsv with at least two columns mapping",
                               "variant IDs between the source callset (first",
                               "column) and the target callset (second column)",
                               "or NA for no match"))
parser$add_argument("--source-site-metrics", metavar=".bed", required=TRUE,
                    help=paste("Site-level summary metrics for all variants",
                               "from the source callset included in --variant-map"))
parser$add_argument("--sample-map", metavar=".tsv", required=TRUE,
                    help=paste("Two-column .tsv mapping sample IDs between the",
                               "source callset (first column) and the target",
                               "callset (second column)"))
parser$add_argument("--out-prefix", metavar="path", type="character",
                    default="./gt_comparison", help="Prefix/path for all output files")
parser$add_argument("--source-gt-dir", metavar="path", default="./",
                    help="Location to search for source sample genotypes [default: pwd]")
parser$add_argument("--target-gt-dir", metavar="path", default="./",
                    help="Location to search for target sample genotypes [default: pwd]")
parser$add_argument("--gt-tsv-suffix", metavar="string", default=".gt.tsv.gz",
                    help=paste("File suffix to assume when searching for sample",
                               "genotype .tsvs; will be appended with sample ID"))
parser$add_argument("--report-by-genotype", action="store_true",
                    help=paste("Split summary results by genotype [default:",
                               "summarize by zygosity]"))
parser$add_argument("--common-af", metavar="float", default=0.01, type="numeric",
                    help="Allele frequency threshold for common variants")
args <- parser$parse_args()

# # DEV:
# args <- list("variant_map" = "~/Downloads/gt_bench_dev_data/revised.vid_map.tsv.gz",
#              "source_site_metrics" = "~/Downloads/gt_bench_dev_data/source.metrics.bed.gz",
#              "sample_map" = "~/Downloads/gt_bench_dev_data/sample.map.tsv",
#              "source_gt_dir" = "~/Downloads/gt_bench_dev_data/",
#              "target_gt_dir" = "~/Downloads/gt_bench_dev_data/",
#              "gt_tsv_suffix" = ".gt.tsv.gz",
#              "report_by_genotype" = TRUE,
#              "common_af" = 0.01,
#              "out_prefix" = "~/scratch/gt_comparison_dev")

# Load input data
vid.map <- load.vid.map(args$variant_map, args$source_site_metrics, args$common_af)

# Load sample map
sid.map <- load.sid.map(args$sample_map)

# Process each sample in parallel
all.sample.res <- lapply(1:nrow(sid.map), function(sidx){
  source.sid <- sid.map[sidx, "source_sid"]
  source.gt.path <- paste(args$source_gt_dir, "/", source.sid,
                          args$gt_tsv_suffix, sep="")
  source.gts <- load.gts(source.gt.path, prefix="source",
                         elig.vids=vid.map$source_vid)
  target.sid <- sid.map[sidx, "target_sid"]
  target.gt.path <- paste(args$target_gt_dir, "/", target.sid,
                          args$gt_tsv_suffix, sep="")
  target.gts <- load.gts(target.gt.path, prefix="target",
                         elig.vids=vid.map$target_vid)
  if(!is.null(source.gts) & !is.null(target.gts)){
    sample.res <- benchmark.gts(source.gts, target.gts, vid.map,
                                args$report_by_genotype)
    cbind(rep(source.sid, nrow(sample.res)), sample.res)
  }else{
    NULL
  }
})
if(length(all.sample.res) > 0){
  res.df <- as.data.frame(do.call("rbind", all.sample.res))
}else{
  res.df <- as.data.frame(matrix(nrow=0, ncol=8))
}
colnames(res.df) <- c("#sample", "class", "subclass", "freq_bin",
                      if(args$report_by_genotype){"genotype"}else{"zygosity"},
                      "no_match", "carrier_match", "gt_match")

# Write results to compressed distribution .tsv
write.table(res.df, paste(args$out_prefix, "gt_comparison.distrib.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

