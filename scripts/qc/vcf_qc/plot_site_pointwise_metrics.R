#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot pointwise variant site-level statistics (as opposed to summary metrics;
# for those, see plot_site_summary_metrics.R)

# For plotting practicality, this script subsets to common variant sites
# if that is not already done by the user upstream


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
require(viridis, quietly=TRUE)
load.constants("all")


##################
# Data Functions #
##################
# Load a site metrics BED file and subset to minimal required columns
read.bed <- function(bed.in, common_af=0, autosomes.only=T){
  df <- read.table(bed.in, header=T, sep="\t", comment.char="",
                   check.names=F, quote="")
  if(autosomes.only){
    df <- df[df[, 1] %in% c(1:22, paste("chr", c(1:22), sep="")), ]
  }
  keep.cols <- c("af", "freq_het", "freq_hom", "hwe")
  df <- as.data.frame(apply(df[, keep.cols], 2, as.numeric))
  df <- df[which(df$af >= common_af), ]
  df <- df[complete.cases(df), ]
  df[, c("hwe.x", "hwe.y")] <- t(apply(df[, c("freq_het", "freq_hom")], 1,
                                       function(freqs){calc.hwe.xy(freqs[1], freqs[2])}))
  return(as.data.frame(df))
}


######################
# Plotting Functions #
######################
# HWE ternary plot
hwe.plot <- function(df, title="All variants", pt.cex=NULL, pt.pch=NULL,
                     pt.alpha=NULL, parmar=c(1.75, 2, 1.25, 1)){
  # Get plot data
  df <- df[which(!is.na(df$hwe.x) & !is.na(df$hwe.y) & !is.na(df$hwe)), ]
  ytop <- sin(60 * pi / 180)
  bonf.p <- 0.05 / nrow(df)
  n.p.bins <- floor(10 * -log10(bonf.p))

  # Dynamically determine point properties
  pt.params <- scatterplot.point.params(nrow(df), cex.start=0.3)
  if(is.null(pt.cex)){
    pt.cex <- pt.params$cex
  }
  if(is.null(pt.pch)){
    pt.pch <- pt.params$pch
  }
  if(is.null(pt.alpha)){
    pt.alpha <- pt.params$alpha
  }
  p.pal <- rev(viridis(n.p.bins + 1, begin=0.1, end=0.95))
  pt.cols <- p.pal[sapply(floor(10 * -log10(df$hwe)),
                          function(v){min(v, n.p.bins)}) + 1]
  pt.cols <- adjustcolor(pt.cols, alpha=pt.alpha)

  # Prep plot area
  prep.plot.area(c(0, 1), c(0, 1), parmar)
  segments(x0=c(0, 0.5, 1), x1=c(0.5, 1, 0),
           y0=c(0, ytop, 0), y1=c(ytop, 0, 0),
           col=annotation.color)

  # Add bottom X axis
  clean.axis(1, at=seq(0, 1, 0.25), label.units="percent", label.line=-1,
             title="Homozygotes", title.line=-0.25)

  # Add left diagonal X axis
  segments(x0=0, x1=0.5, y0=0, y1=ytop)
  diag.ax.spacing  <- seq(0, 1, 0.25)
  diag.ax.x.at <- diag.ax.spacing * cos(60 * pi / 180)
  diag.ax.y.at <- diag.ax.spacing * sin(60 * pi / 180)
  diag.ax.tck <- 0.025
  diag.ax.x.adj <- diag.ax.tck * cos(30 * pi / 180)
  diag.ax.y.adj <- diag.ax.tck * sin(30 * pi / 180)
  segments(x0=diag.ax.x.at, x1=diag.ax.x.at - diag.ax.x.adj,
           y0=diag.ax.y.at, y1=diag.ax.y.at + diag.ax.y.adj, xpd=T)
  text(x=diag.ax.x.at+0.025, y=diag.ax.y.at+0.015, pos=2, cex=4.5/6,
       labels=paste(seq(0, 100, 25), "%", sep=""), xpd=T)
  text(x=mean(diag.ax.x.at) - (10*diag.ax.x.adj),
       y=mean(diag.ax.y.at) + (10*diag.ax.y.adj),
       srt=60, labels="Heterozygotes", xpd=T)

  # Add legend
  legend.left <- 0.75
  legend.right <- 0.8
  legend.y.breaks <- seq(2*ytop/3, ytop, length.out=4)
  legend.col.breaks <- c(0, floor(10 * -log10(0.05)),
                         floor(10 * -log10(bonf.p)), length(p.pal))
  hwe.leg.labs <- list(format.pval(0.05, equality=">="),
                       format.pval(0.05, equality="<"),
                       format.pval(bonf.p, equality="<"))
  text(x=legend.left+0.025, y=mean(legend.y.breaks[3:4]),
       pos=2, labels=bquote(italic(P)))
  sapply(1:3, function(i){
    k.range <- legend.col.breaks[i]:legend.col.breaks[i+1]
    k.rect.y.breaks <- seq(legend.y.breaks[i], legend.y.breaks[i+1],
                           length.out=length(k.range))
    rect(xleft=rep(legend.left, length(k.range)),
         xright=rep(legend.right, length(k.range)),
         ybottom=k.rect.y.breaks[-length(k.rect.y.breaks)],
         ytop=k.rect.y.breaks[-c(1)], col=p.pal[k.range],
         border=p.pal[k.range], lwd=0.5)
    text(x=legend.left-0.015,
         y=mean(legend.y.breaks[c(i, i+1)]),
         pos=4, cex=5/6, xpd=T,
         labels=if(i==1){
           expression("" > 0.05)
         }else if(i==2){
           expression("" <= 0.05)
         }else if(i==3){
           expression("" < "Bonf.")
         })
  })
  rect(xleft=rep(legend.left, 3), xright=rep(legend.right, 3),
       ybottom=legend.y.breaks[1:3], ytop=legend.y.breaks[2:4],
       col=NA, xpd=T, border="white")
  rect(xleft=legend.left, xright=legend.right,
       ybottom=legend.y.breaks[1], ytop=legend.y.breaks[4],
       col=NA)

  # Add title & subtitle
  mtext(3, line=0.4, text=title)
  n.all <- nrow(df)
  n.pass <- sum(df$hwe >= bonf.p)
  n.formatted <- clean.numeric.labels(c(n.pass, n.all), acceptable.decimals=2)
  n.formatted[1] <- gsub("k|M|B|T", "", n.formatted[1])
  pct.pass <- n.pass / n.all
  subtitle <- paste(n.formatted[1], " / ",
                    n.formatted[2], " (", round(100*pct.pass, 1),
                    "%) pass HWE", sep="")
  mtext(3, cex=5/6, text=subtitle, line=-0.4)

  # Add points
  points(df$hwe.x, df$hwe.y, pch=1, cex=pt.cex, col=pt.cols, xpd=T)

  return(c(pct.pass, n.all))
}

# Plot wrapper function for all site-level pointwise plots
pointwise.plots <- function(df, out.prefix, fname.suffix="all",
                            title="All variants"){

  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.prefix <- gsub("[ ]+", "_", fname.suffix)

  # HWE plot as .png
  png(paste(out.prefix, fname.suffix, "hwe.png", sep="."),
      height=2.25*300, width=2.25*300, res=300)
  m.tmp <- hwe.plot(df, title=title)
  dev.off()

  ss.df[1, ] <- c(paste(ss.prefix, "common_hwe", sep="."), "pct_pass", m.tmp)

  return(ss.df)
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot pointwise variant metrics")
parser$add_argument("--snvs", metavar=".bed", type="character",
                    help="SNV site metrics from clean_site_metrics.py")
parser$add_argument("--indels", metavar=".bed", type="character",
                    help="Indel site metrics from clean_site_metrics.py")
parser$add_argument("--svs", metavar=".bed", type="character",
                    help="SV site metrics from clean_site_metrics.py")
parser$add_argument("--combine", action="store_true", default=FALSE,
                    help="Also generate a combined set of plots for all variant types")
parser$add_argument("--common-af", metavar="float", default=0.01, type="numeric",
                    help="Allele frequency threshold for common variants")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV:
# args <- list("snvs" = NULL,
#              "indels" = NULL,
#              "svs" = "~/scratch/YL.sv.site_metrics.dev.sv.sites.common.bed.gz",
#              "combine" = FALSE,
#              "common_af" = 0.01,
#              "out_prefix" = "~/scratch/qc.test")

# Ensure at least one of --snvs, --indels, or --svs is present
if(is.null(args$snvs) & is.null(args$indels) & is.null(args$svs)){
  stop("At least one of --snvs, --indels, or --svs must be provided")
}

# Load & plot SNVs, if provided
if(!is.null(args$snvs)){
  # Load SNV data
  snv.df <- read.bed(args$snvs, args$common_af)

  # Plot SNV metrics
  snv.ss <- pointwise.plots(snv.df, args$out_prefix, fname.suffix="snv",
                  title="Common SNVs")
}else{
  snv.df <- NULL
  snv.ss <- NULL
}

# Load & plot indels, if provided
if(!is.null(args$indels)){
  # Load indel data
  indel.df <- read.bed(args$indels, args$common_af)

  # Plot indel metrics
  indel.ss <- pointwise.plots(indel.df, args$out_prefix, fname.suffix="indel",
                  title="Common indels")

}else{
  indel.df <- NULL
  indel.ss <- NULL
}

# Load & plot SVs, if provided
if(!is.null(args$svs)){
  # Load SV data
  sv.df <- read.bed(args$svs, args$common_af)

  # Plot SV metrics
  sv.ss <- pointwise.plots(sv.df, args$out_prefix, fname.suffix="sv",
                  title="Common SVs")
}else{
  sv.df <- NULL
  sv.ss <- NULL
}

# Combine & plot all variant types, if optioned
if(args$combine & sum(sapply(list(snv.df, indel.df, sv.df), is.null)) < 2){
  # Merge all data
  all.df <- do.call("rbind", list(snv.df, indel.df, sv.df))

  # Plot all metrics
  all.ss <- pointwise.plots(all.df, args$out_prefix, fname.suffix="all",
                  title="All common variants")
}else{
  all.ss <- NULL
}

# Combine summary statistics and write to .tsv
ss.out <- do.call("rbind", list(snv.ss, indel.ss, sv.ss, all.ss))
colnames(ss.out)[1] <- paste("#", colnames(ss.out)[1], sep="")
write.table(ss.out, paste(args$out_prefix, "summary_metrics.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)

