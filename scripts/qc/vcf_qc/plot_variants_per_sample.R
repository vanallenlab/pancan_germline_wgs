#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot distribution of variants per genome


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
# Load and process genotype counts per sample
load.gt.counts <- function(gt.dist.in){
  # Read data
  gt <- read.table(gt.dist.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(gt) <- gsub("#", "", colnames(gt))

  # For purposes of vizualization, lump "other" GTs in with homozygotes
  gt$hom <- gt$hom + gt$other
  gt$other <- NULL

  # Add collapsed parent categories per variant class
  vcs <- unique(gt$class)
  sids  <- unique(gt$sample)
  add.rows <- as.data.frame(do.call("rbind", lapply(vcs, function(vc){
    do.call("rbind", lapply(sids, function(sid){
      sub.df <- gt[which(gt$sample == sid & gt$class == vc), ]
      c(sid, vc, NA, NA, apply(sub.df[, c("het", "hom")], 2, sum))
    }))
  })))
  colnames(add.rows) <- colnames(gt)
  gt <- as.data.frame(rbind(gt, add.rows))
  gt[, c("het", "hom")] <- apply(gt[, c("het", "hom")], 2, as.numeric)

  # Add parent row for total variation per genome across all classes
  total.rows <- as.data.frame(do.call("rbind", lapply(sids, function(sid){
    sub.df <- gt[which(gt$sample == sid & gt$class %in% vcs & is.na(gt$subclass)), ]
    c(sid, "all", NA, NA, apply(sub.df[, c("het", "hom")], 2, sum))
  })))
  colnames(total.rows) <- colnames(gt)
  gt <- as.data.frame(rbind(gt, total.rows))
  gt[, c("het", "hom")] <- apply(gt[, c("het", "hom")], 2, as.numeric)

  return(gt)
}

# Load sample labels and subset to a list of eligible samples
load.labels <- function(tsv.in, elig.samples=NULL){
  if(is.null(tsv.in)){
    return(NULL)
  }
  labs <- read.table(tsv.in, header=F, sep="\t")
  if(!is.null(elig.samples)) {
    labs <- labs[which(labs[, 1] %in% elig.samples), ]
  }
  l.n <- as.character(labs[, 1])
  l.v <- as.character(labs[, 2])
  names(l.v) <- l.n

  return(l.v)
}


######################
# Plotting Functions #
######################
# Scatterplot of sample heterozygosity between two classes of variants
plot.heterozygosity <- function(gt.counts, vc1, vc2, pop=NULL,
                                title="Heterozygosity", stat.lab.hex=0.04,
                                pt.cex=0.65, parmar=c(2, 2.5, 1, 3.5)){
  # Compute heterozygosity for all samples for each variant class
  samples <- unique(gt.counts$sample)
  h.df <- as.data.frame(do.call("cbind", lapply(c(vc1, vc2), function(vc){
    sapply(samples, function(sid){
      sv1 <- gt.counts[(which(gt.counts$sample == sid
                              & gt.counts$class == vc
                              & is.na(gt.counts$subclass))),
                       c("het", "hom")]
      if(length(sv1) == 0){NA}else{sv1$het / sum(sv1)}
    })
  })))
  rownames(h.df) <- samples
  h.df <- h.df[which(complete.cases(h.df)), ]
  colnames(h.df) <- c(vc1, vc2)

  # Reorder samples from smallest to largest pop, if optioned
  if(!is.null(pop)){
    pop.k <- table(pop[rownames(h.df)])
    pop.idx <- length(pop.k):1
    names(pop.idx) <- names(pop.k)[order(pop.k)]
    h.df <- h.df[order(pop.idx[pop[rownames(h.df)]], h.df[, vc1], h.df[, vc2]), ]
  }

  # Get pointwise plotting parameters
  pw.params <- scatterplot.point.params(nrow(h.df))
  if(is.null(pop)){
    pw.col <- rep("black", nrow(h.df))
  }else{
    n.pops <- length(unique(pop))
    pop <- toupper(pop)
    if(all(pop %in% names(pop.colors))){
      pop.pal <- pop.colors[sort(unique(pop))]
    }else{
      pop.pal <- categorical.rainbow(n.pops)
      names() <- unique(pop)
    }
    pw.col <- pop.colors[pop[rownames(h.df)]]
  }

  # Prepare plotting area
  xlims <- ylims <- range(h.df[, c(vc1, vc2)])
  prep.plot.area(xlims, ylims, parmar, xaxs="r", yaxs="r")
  abline(0, 1, col=annotation.color)

  # Add points
  points(h.df, pch=pw.params$pch, cex=0.75*pw.params$cex, xpd=T,
         col=sapply(pw.col, adjustcolor, alpha=pw.params$alpha))

  # Add axes & title
  clean.axis(1, label.units="percent", infinite=T, max.ticks=5,
             title=paste(var.class.abbrevs[vc1], "s", sep=""),
             label.line=-0.85, title.line=0)
  clean.axis(2, label.units="percent", infinite=T, max.ticks=5,
             title=paste(var.class.abbrevs[vc2], "s", sep=""),
             label.line=-0.75, title.line=0.5)
  mtext(3, text=title)

  # Add legends
  r2 <- cor(h.df[, vc1], h.df[, vc2])^2
  text(x=par("usr")[1]-(stat.lab.hex*diff(par("usr")[1:2])),
       y=par("usr")[4]-(2.5*stat.lab.hex*diff(par("usr")[3:4])),
       pos=4, cex=5/6,
       labels=bquote(italic(R)^2*"="*.(formatC(round(r2, 2), digits=2))))
  if(!is.null(pop)){
    pop.labs.at <- seq(par("usr")[1], par("usr")[2],
                       length.out=2+n.pops)[-c(1, n.pops + 2)]
    if(all(pop %in% names(pop.abbreviations))){
      pop.labs <- pop.abbreviations[names(pop.pal)]
    }else{
      pop.labs <- names(pop.pal)
    }
    yaxis.legend(pop.labs, x=par("usr")[2], y.positions=rev(pop.labs.at),
                 sep.wex=0.05*diff(par("usr")[1:2]), colors=pop.pal)
  }

  return(c(r2, nrow(h.df)))
}

# Wrapper function to handle inter-class comparison plots
plot.vc.comparisons <- function(gt.counts, vc1, vc2, out.prefix, pop=NULL){

  ss.df <- data.frame("analysis"=character(), "measure"=character(),
                      "value"=numeric(), "n"=numeric())
  ss.prefix <- paste(vc1, "vs", vc2, sep="_")

  # Scatterplot of heterozygosity per sample between all pairs of variant classes
  png(paste(out.prefix, ss.prefix, "heterozygosity.png", sep="."),
      height=2.25*300, width=2.75*300, res=300)
  m.tmp <- plot.heterozygosity(gt.counts, vc1, vc2, pop=pop)
  dev.off()

  ss.df[1, ] <- c(paste(ss.prefix, "heterozygosity_cor", sep="."), "r2", m.tmp)

  # Scatterplot of total variant counts per sample between all pairs of variant classes
  # TODO: implement this

  return(ss.df)
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot variant distributions per genome")
parser$add_argument("--genotype-dist-tsv", required=TRUE, metavar=".tsv",
                    help=paste("Input .tsv of genotype counts per sample ",
                               "generated by clean_sample_genotypes.py"))
parser$add_argument("--ancestry-labels", metavar=".tsv",
                    help=paste("Optional two-column .tsv mapping sample IDs to ",
                               "ancestry labels. If provided, all plots will ",
                               "distinguish between ancestries."))
parser$add_argument("--phenotype-labels", metavar=".tsv",
                    help=paste("Optional two-column .tsv mapping sample IDs to ",
                               "phenotype labels. If provided, some plots will ",
                               "group samples by phenotype for easier comparisons."))
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV:
# args <- list("genotype_dist_tsv" = "~/scratch/gt_plot_dev/dfci-g2c.v1.initial_qc.chr22.genotype_distribution.merged.tsv.gz",
#              "ancestry_labels" = "~/scratch/dfci-g2c.v1.qc_ancestry.tsv",
#              "phenotype_labels" = "~/scratch/dfci-g2c.v1.qc_phenotype.tsv",
#              "out_prefix" = "~/scratch/g2c_vcf_qc_dev")

# Load sample genotype counts
gt.counts <- load.gt.counts(args$genotype_dist_tsv)
samples <- unique(gt.counts$sample)

# Load sample ancestry and phenotypes, if optioned
pop <- load.labels(args$ancestry_labels, samples)
pheno <- load.labels(args$ancestry_labels, samples)

# Triple waterfall plot of counts per sample by class
# TODO: implement this

# Inter-class comparisions across samples
ic.ss <- lapply(list(c("snv", "indel"), c("snv", "sv"), c("indel", "sv")), function(vcs){
  vc1 <- vcs[1]; vc2 <- vcs[2]
  if(all(vcs %in% gt.counts$class)){
    plot.vc.comparisons(gt.counts, vc1, vc2, args$out_prefix, pop)
  }
})
if(length(ic.ss) > 0){
  ic.ss <- as.data.frame(do.call("rbind", ic.ss))
}
