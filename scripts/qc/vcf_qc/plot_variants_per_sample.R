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
# Helper function to plot a single count waterfall; called within main.waterfall()
plot.count.waterfall <- function(gt.counts, vc, pop=NULL, pheno=NULL,
                                 samples.sorted=NULL, pop.space.wex=0.02,
                                 pheno.space.wex=0.005, parmar=c(0.25, 2.75, 0.5, 0.15)){
  # Format plotting data
  counts <- gt.counts[which(gt.counts$class == vc
                            & is.na(gt.counts$subclass)),
                      c("sample", "hom", "het")]
  rownames(counts) <- counts$sample
  if(!is.null(samples.sorted)){
    counts <- counts[samples.sorted, ]
  }
  counts$sample <- NULL
  s.feats <- data.frame("order" = 1:nrow(counts), row.names=rownames(counts))
  s.feats$pop <- if(is.null(pop)){NA}else{pop[rownames(counts)]}
  s.feats$pheno <- if(is.null(pop)){NA}else{pheno[rownames(counts)]}
  pop.breaks.x <- which(sapply(2:nrow(s.feats), function(r){
    s.feats$pop[r] != s.feats$pop[r-1]
  }))
  pheno.breaks.x <- which(sapply(2:nrow(s.feats), function(r){
    (s.feats$pop[r] == s.feats$pop[r-1]
     & s.feats$pheno[r] != s.feats$pheno[r-1])
  }))
  n.samples <- length(rownames(counts))

  # Get plot parameters
  xlims <- c(-0.01*n.samples,
             n.samples * (1 + (pop.space.wex * (length(pop.breaks.x) - 1))
                          + (pheno.space.wex * (length(pheno.breaks.x) - 1))))
  ylims <- c(0, max(apply(counts, 1, sum, na.rm=T)))
  hom.col <- adjust.color.hsb(var.class.colors[vc], s=-0.01, b=-0.01)
  het.col <- adjust.color.hsb(var.class.colors[vc], s=0.025, b=0.075)

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar)

  # Add bars
  prev.x <- 0
  for(i in 1:n.samples){
    if(i %in% pop.breaks.x){
      prev.x <- prev.x + (pop.space.wex * n.samples)
    }else if(i %in% pheno.breaks.x){
      prev.x <- prev.x + (pheno.space.wex* n.samples)
    }
    rect(xleft=prev.x, xright=prev.x+1, ybottom=c(0, counts[i, 1]),
         ytop=cumsum(as.numeric(counts[i, 1:2])),
         col=c(hom.col, het.col), border=c(hom.col, het.col), lwd=0.2)
    prev.x <- prev.x + 1
  }

  # Add axis
  clean.axis(2, title=paste(var.class.abbrevs[vc], "s", sep=""),
             label.units="count", infinite.positive=T, max.ticks=4,
             title.line=0.75)
}

# Wrapper for main waterfall plot
main.waterfall <- function(gt.counts, out.prefix, pop=NULL, pheno=NULL,
                           pop.space.wex=0.03, pheno.space.wex=0.005){
  # Determine number of panels and figure sizing
  vcs <- intersect(names(var.class.names), unique(gt.counts$class))
  n.panels <- length(vcs)
  has.pop <- !is.null(pop)
  has.pheno <- !is.null(pheno)
  bottom.tracks <- sum(c(has.pop, has.pheno))
  samples <- unique(gt.counts$sample)

  # Compute each sample's total variant counts
  total.counts <- as.data.frame(do.call("cbind", lapply(vcs, function(vc){
    sapply(samples, function(sid){
      s.v <- as.numeric(gt.counts[which(gt.counts$sample == sid
                                        & gt.counts$class == vc
                                        & is.na(gt.counts$subclass)), c("het", "hom")])
      if(length(s.v) == 0){0}else{sum(s.v, na.rm=T)}
    })
  })))
  colnames(total.counts) <- vcs
  rownames(total.counts) <- samples

  # Determine sample order as follows:
  # 1. Group samples on X axis by population from most to least total variants
  # 2. Order samples within each group first by phenotype (largest to smallest N)
  # 3. Order samples within `2` from most to least variation
  if(has.pop){
    pop <- pop[samples]
    pops <- unique(pop)
    pop.m <- sapply(pops, function(p){
      median(apply(total.counts[names(pop[which(pop == p)]), ], 1, sum))
    })
    pop.rank <- 1:length(pop.m)
    names(pop.rank) <- names(pop.m)[order(pop.m, decreasing=TRUE)]
  }
  if(has.pheno){
    pheno <- pheno[samples]
    pheno.k <- table(pheno)
    pheno.rank <- 1:length(pheno.k)
    names(pheno.rank) <- names(pheno.k)[order(pheno.k, decreasing=TRUE)]
  }
  s.z <- apply(apply(total.counts, 2, scale), 1, mean)
  s.rank <- 1:length(samples)
  names(s.rank) <- rownames(total.counts)[order(s.z, decreasing=T)]
  order.df <- data.frame("counts" = s.rank, row.names=names(s.rank))
  if(has.pop){
    order.df$pop <- pop.rank[pop[rownames(order.df)]]
  }
  if(has.pheno){
    order.df$pheno <- pheno.rank[pheno[rownames(order.df)]]
  }
  order.df <- order.df[, intersect(c("pop", "pheno", "counts"),
                                   colnames(order.df))]
  samples.sorted <- rownames(order.df)[do.call(order, order.df)]

  # Define plot parameters
  # Target dimensions for 3 rows + 2 marker fields: 6.8" wide x 4.2" tall
  pdf.width <- 6.8
  pdf.height <- (3*0.42*n.panels) + (bottom.tracks*0.42/2)

  # Generate waterfall plot
  pdf(paste(out.prefix, "variants_per_genome.waterfall.pdf", sep="."),
      height=pdf.height, width=pdf.width)
  layout(matrix(1:(n.panels+(bottom.tracks > 0)), byrow=T, ncol=1),
         heights=c(rep(3, n.panels), if(bottom.tracks > 0){1}))
  sapply(vcs, function(vc){
    plot.count.waterfall(gt.counts, vc, pop, pheno, samples.sorted,
                         pop.space.wex, pheno.space.wex)
  })
  dev.off()

}

# Helper function to handle scatterplots used for inter-class comparisons
interclass.scatter <- function(plot.df, pop=NULL, title=NULL, label.units=NULL,
                               xlims=NULL, ylims=NULL, diag.lm=F,
                               stat.lab.hex=0.04, pt.cex=0.65,
                               parmar=c(2, 2.5, 1, 3.5)){
  # Get pointwise plotting parameters
  pw.params <- scatterplot.point.params(nrow(plot.df))
  if(is.null(pop)){
    pw.col <- rep("black", nrow(plot.df))
  }else{
    n.pops <- length(unique(pop))
    pop <- toupper(pop)
    if(all(pop %in% names(pop.colors))){
      pop.pal <- pop.colors[sort(unique(pop))]
    }else{
      pop.pal <- categorical.rainbow(n.pops)
      names() <- unique(pop)
    }
    pw.col <- pop.colors[pop[rownames(plot.df)]]
  }

  # Prepare plotting area
  if(is.null(xlims)){
    xlims <- range(plot.df[, 1])
  }
  if(is.null(ylims)){
    ylims <- range(plot.df[, 2])
  }
  prep.plot.area(xlims, ylims, parmar, xaxs="r", yaxs="r")
  if(diag.lm){abline(0, 1, col=annotation.color)}

  # Add points
  points(plot.df, pch=pw.params$pch, cex=0.75*pw.params$cex, xpd=T,
         col=sapply(pw.col, adjustcolor, alpha=pw.params$alpha))

  # Add axes & title
  mtext(3, text=title)
  clean.axis(1, label.units=label.units, infinite=T, max.ticks=5,
             title=tryCatch(paste(var.class.abbrevs[colnames(plot.df)[1]], "s", sep=""),
                            error=function(e){"X"}),
             label.line=-0.85, title.line=0)
  clean.axis(2, label.units=label.units, infinite=T, max.ticks=5,
             title=tryCatch(paste(var.class.abbrevs[colnames(plot.df)[2]], "s", sep=""),
                            error=function(e){"Y"}),
             label.line=-0.75, title.line=0.5)

  # Add title & legends
  r2 <- cor(plot.df[, 1], plot.df[, 2])^2
  text(x=par("usr")[1]-(stat.lab.hex*diff(par("usr")[1:2])),
       y=par("usr")[4]-(2.5*stat.lab.hex*diff(par("usr")[3:4])),
       pos=4, cex=5/6,
       labels=bquote(italic(R)^2*"="*.(formatC(round(r2, 2), digits=2))))
  if(!is.null(pop)){
    pop.labs.at <- seq(par("usr")[3], par("usr")[4],
                       length.out=2+n.pops)[-c(1, n.pops + 2)]
    if(all(pop %in% names(pop.abbreviations))){
      pop.labs <- pop.abbreviations[names(pop.pal)]
    }else{
      pop.labs <- names(pop.pal)
    }
    yaxis.legend(pop.labs, x=par("usr")[2], y.positions=rev(pop.labs.at),
                 sep.wex=0.05*diff(par("usr")[1:2]), colors=pop.pal, lwd=4)
  }

  return(r2)
}

# Scatterplot of sample heterozygosity between two classes of variants
plot.heterozygosity <- function(gt.counts, vc1, vc2, pop=NULL,
                                title="Heterozygosity"){
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

  # Generate plot
  xlims <- ylims <- range(h.df[, c(vc1, vc2)])
  r2 <- interclass.scatter(h.df, pop, title, label.units="percent",
                           xlims=xlims, ylims=ylims, diag.lm=T)

  return(c(r2, nrow(h.df)))
}

# Scatterplot of overall variant counts between two classes of variants
plot.count.comparisons <- function(gt.counts, vc1, vc2, pop=NULL,
                                   title="Variant count"){
  # Compute total number of variants for all samples for each variant class
  samples <- unique(gt.counts$sample)
  c.df <- as.data.frame(do.call("cbind", lapply(c(vc1, vc2), function(vc){
    sapply(samples, function(sid){
      sum(gt.counts[(which(gt.counts$sample == sid
                           & gt.counts$class == vc
                           & is.na(gt.counts$subclass))),
                    c("het", "hom")])
    })
  })))
  rownames(c.df) <- samples
  c.df <- c.df[which(complete.cases(c.df)), ]
  colnames(c.df) <- c(vc1, vc2)

  # Reorder samples from smallest to largest pop, if optioned
  if(!is.null(pop)){
    pop.k <- table(pop[rownames(c.df)])
    pop.idx <- length(pop.k):1
    names(pop.idx) <- names(pop.k)[order(pop.k)]
    c.df <- c.df[order(pop.idx[pop[rownames(c.df)]], c.df[, vc1], c.df[, vc2]), ]
  }

  # Generate plot
  r2 <- interclass.scatter(c.df, pop, title, label.units="count")

  return(c(r2, nrow(c.df)))
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
  png(paste(out.prefix, ss.prefix, "variant_counts.png", sep="."),
      height=2.25*300, width=2.75*300, res=300)
  m.tmp <- plot.count.comparisons(gt.counts, vc1, vc2, pop=pop)
  dev.off()

  ss.df[2, ] <- c(paste(ss.prefix, "variant_count_cor", sep="."), "r2", m.tmp)

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
pheno <- load.labels(args$phenotype_labels, samples)

# Triple waterfall plot of counts per sample by class
main.waterfall(gt.counts, args$out_prefix, pop, pheno)

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
