#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot variant site-level summary distributions for quality control


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2CR, quietly=TRUE)
require(viridis, quietly=TRUE)
load.constants("all")

# Declare local custom constants
sv.dens.bw <- c("DEL" = 0.5,
                "DUP" = 4/5,
                "CNV" = 0.5,
                "INS" = 10,
                "INV" = 3/4,
                "CPX" = 1/3)


##################
# Data Functions #
##################
# Clean size or AF breaks encoded as strings
clean.breaks <- function(vals){
  as.numeric(gsub("ge|lt|bp", "", vals))
}

# Import a precomputed compressed summary distribution
read.distrib <- function(tsv.in, key.cols=1:2){
  if(is.null(tsv.in)){return(NULL)}
  df <- read.table(tsv.in, header=T, sep="\t", comment.char="",
                   quote="", check.names=F)
  colnames(df)[1] <- gsub("#", "", colnames(df)[1])
  df[, -key.cols] <- apply(df[, -key.cols], 2, as.numeric)
  dist.breaks <- clean.breaks(colnames(df)[-key.cols])
  list("df" = df, "breaks" = dist.breaks)
}

# Load all SV sizes from an SV sites BED and split by subclass
read.sv.sizes <- function(tsv.in){
  if(is.null(tsv.in)){return(NULL)}
  df <- read.table(tsv.in, header=T, sep="\t", quote="",
                   check.names=F, comment.char="")[, c("subclass", "size")]
  sv.classes.in.sites <- setdiff(intersect(names(sv.colors), unique(df$subclass)),
                                 c("CTX", "BND"))
  sv.sizes <- lapply(sv.classes.in.sites, function(vsc){
    log10(as.numeric(df$size[which(df$subclass == vsc)]))
  })
  names(sv.sizes) <- sv.classes.in.sites
  return(sv.sizes)
}


######################
# Plotting Functions #
######################
# Summary plot of variant counts by class & subclass
plot.counts.by.vsc <- function(df, bar.sep=0.1, parmar=c(0.1, 7, 2, 2)){
  # Simplify count data
  k <- log10(apply(df[, -c(1:2)], 1, sum))
  k.order <- order(k)
  k <- k[k.order]
  df <- df[k.order, ]
  names(k) <- df$subclass

  # Get plot parameters
  xlims <- c(0, max(ceiling(k)))
  ylims <- c(0, length(k)+bar.sep)
  bar.cols <- var.class.colors[df$class]
  bar.labs <- sapply(10^k, clean.numeric.labels)

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar)
  # abline(v=0, col="gray90", lwd=2)

  # Add bars & count labels
  rect(xleft=rep(0, length(k)), xright=k,
       ybottom=(1:length(k)) - 1 + bar.sep, ytop=(1:length(k)) - bar.sep,
       col=bar.cols, xpd=T)
  text(x=k-0.5, y=(1:length(k))-0.55, pos=4, cex=5/6, labels=bar.labs, xpd=T)

  # Add top Y axis
  minor.ticks <- log10(logscale.minor)
  minor.ticks <- minor.ticks[which(minor.ticks < max(ylims))]
  axis(3, at=minor.ticks, labels=NA, tck=-0.01, lwd=0.75)
  major.at <- log10(logscale.major.bp)
  major.labels <- sapply(logscale.major.bp, clean.numeric.labels)
  major.labels[which(!major.labels %in% c("1", "1k", "1M", "1B"))] <- NA
  clean.axis(3, at=major.at, labels=major.labels, title="Total variants",
             label.line=-0.8, title.line=-0.1, tck=-0.03,
             infinite.positive=TRUE)

  # Add left margin labels
  axis(2, at=(1:length(k))-0.5, las=2, line=-0.9, tick=F,
       labels=var.subclass.names.short[df$subclass], cex.axis=5/6)
  indel.vc.x <- sv.vc.x <- -5.75
  indel.bracket.y <- c(min(which(df$class == "indel")) - 1 + bar.sep,
                       max(which(df$class == "indel")) - bar.sep)
  staple.bracket(x0=indel.vc.x, x1=indel.vc.x,
                 y0=indel.bracket.y[1], y1=indel.bracket.y[2])
  text(x=indel.vc.x+0.5, y=mean(indel.bracket.y)-0.1, labels="Indels\n(1-49 bp)",
       cex=5/6, pos=2, xpd=T)

  sv.bracket.y <- c(min(which(df$class == "sv")) - 1 + bar.sep,
                    max(which(df$class == "sv")) - bar.sep)
  staple.bracket(x0=sv.vc.x, x1=sv.vc.x,
                 y0=sv.bracket.y[1], y1=sv.bracket.y[2])
  text(x=sv.vc.x+0.5, y=mean(sv.bracket.y)-0.1,
       labels="Structural\nvariants\n(>49 bp)",
       cex=5/6, pos=2, xpd=T)
}

# Volcano of signed variant sizes
plot.size.volcano <- function(size.d, snv.width=0.2, snv.gap=0.15, indel.gap=0,
                              minor.tck=-0.01, major.tck=-0.0275,
                              parmar=c(2, 2.6, 0.25, 0.1)){
  # Get partitioned densities
  df <- size.d$df
  breaks <- log10(size.d$breaks)
  snv.k <- log10(sum(as.numeric(df[which(df$class == "snv"), -c(1:2)])))
  ins.k <- log10(as.numeric(df[which(df$subclass == "ins"), -c(1:2)]))
  del.k <- log10(as.numeric(df[which(df$subclass == "del"), -c(1:2)]))
  GAIN.k <- log10(as.numeric(apply(df[which(df$subclass %in% c("INS", "DUP", "CNV")), -c(1:2)], 2, sum)))
  LOSS.k <- log10(as.numeric(df[which(df$subclass == "DEL"), -c(1:2)]))

  # Get plot values
  xlims <- c(-1, 1) * (max(breaks) + 0.1 + (0.5*snv.width) + snv.gap + indel.gap)
  xlims[1] <- xlims[1]-0.1
  ylims <- c(0, log10(max(df[, -c(1:2)])))

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="r")

  # Add SNVs
  rect(xleft=-0.5*snv.width, xright=0.5*snv.width,
       ybottom=0, ytop=snv.k, col=var.class.colors["snv"],
       border=NA, bty="n")
  # segments(x0=-0.5*class.gap, x1=0.5*class.gap,
  #          y0=snv.k, y1=snv.k, xpd=T)
  # axis(1, at=c(-0.5, 0.5) * snv.width, tck=0, labels=NA)
  # axis(1, at=0, tick=F, line=-0.9, cex.axis=5/6, labels=0)

  # Add indel polygons
  indel.idx <- which(breaks >= 0 & breaks < log10(50))
  ins.xy <- step.function(x=breaks[indel.idx] + snv.gap + (0.5 * snv.width),
                          y=ins.k[indel.idx], offset=1)
  polygon(x=c(ins.xy$x, rev(ins.xy$x)),
          y=c(ins.xy$y, rep(0, length(ins.xy$x))),
          col=var.class.colors["indel"], border=NA, bty="n")
  # points(x=ins.xy$x, y=ins.xy$y, type="l")
  del.xy <- step.function(x=breaks[indel.idx] + snv.gap + (0.5 * snv.width),
                          y=del.k[indel.idx], offset=1)
  polygon(x=-c(del.xy$x, rev(del.xy$x)),
          y=c(del.xy$y, rep(0, length(del.xy$x))),
          col=var.class.colors["indel"], border=NA, bty="n")
  # points(x=-del.xy$x, y=del.xy$y, type="l")


  # Add indel axes
  logscale.indels <- logscale.minor[which(logscale.minor >= 1 & logscale.minor < 50)]
  axis(1, at=log10(logscale.indels) + snv.gap + (0.5 * snv.width),
       tck=minor.tck, labels=NA, lwd=0.75)
  axis(1, at=log10(c(1, 50)) + snv.gap + (0.5 * snv.width), labels=NA, tck=0)
  axis(1, at=-log10(logscale.indels) - snv.gap - (0.5 * snv.width),
       tck=minor.tck, labels=NA, lwd=0.75)
  axis(1, at=-log10(c(1, 50)) - snv.gap - (0.5 * snv.width), labels=NA, tck=0)
  clean.axis(1, at=0:1 + snv.gap + (0.5 * snv.width), tck=major.tck,
             labels=c(1, 10), cex.axis=4.5/6, label.line=-0.9)
  clean.axis(1, at=0:-1 - snv.gap - (0.5 * snv.width), tck=major.tck,
             labels=c(1, 10), cex.axis=4.5/6, label.line=-0.9)

  # Add SV polygons
  sv.idx <- which(breaks >= log10(50))
  gain.xy <- step.function(x=breaks[sv.idx] + snv.gap + indel.gap + (0.5*snv.width),
                           y=GAIN.k[sv.idx], offset=1)
  polygon(x=c(gain.xy$x, rev(gain.xy$x)),
          y=c(gain.xy$y, rep(0, length(gain.xy$x))),
          col=var.class.colors["sv"], border=NA, bty="n", xpd=T)
  # points(x=gain.xy$x, y=gain.xy$y, type="l", xpd=T)
  loss.xy <- step.function(x=breaks[sv.idx] + snv.gap + indel.gap + (0.5*snv.width),
                           y=LOSS.k[sv.idx], offset=1)
  polygon(x=-c(loss.xy$x, rev(loss.xy$x)),
          y=c(loss.xy$y, rep(0, length(loss.xy$x))),
          col=var.class.colors["sv"], border=NA, bty="n", xpd=T)
  # points(x=-loss.xy$x, y=loss.xy$y, type="l", xpd=T)

  # Add SV axes
  # logscale.major.bp.labels[which(logscale.major.bp.labels == "1Mb")] <- ">1Mb"
  # logscale.major.bp.labels[which(logscale.major.bp.labels == "100kb")] <- NA
  logscale.major.bp.labels <- gsub("b|p", "", logscale.major.bp.labels)
  logscale.svs <- logscale.minor[which(logscale.minor >= 50 & logscale.minor <= 1e6)]
  sv.lbp.idx <- which(logscale.major.bp >= 50 & logscale.major.bp <= 1e6)
  sv.lab.at <- logscale.major.bp[sv.lbp.idx]
  axis(1, at=log10(logscale.svs) + snv.gap + indel.gap + (0.5*snv.width),
       tck=minor.tck, labels=NA, lwd=0.75)
  axis(1, at=-log10(logscale.svs) - snv.gap - indel.gap - (0.5*snv.width),
       tck=minor.tck, labels=NA, lwd=0.75)
  clean.axis(1, at=log10(logscale.major.bp[sv.lbp.idx]) + snv.gap + indel.gap + (0.5*snv.width),
             tck=major.tck, cex.axis=4.5/6, label.line=-0.9,
             labels=logscale.major.bp.labels[sv.lbp.idx])
  clean.axis(1, at=-log10(logscale.major.bp[sv.lbp.idx]) - snv.gap - indel.gap - (0.5*snv.width),
             tck=major.tck, cex.axis=4.5/6, label.line=-0.9,
             labels=logscale.major.bp.labels[sv.lbp.idx])

  # Add directional X axis titles
  axis(1, at=mean(c(-snv.gap - (0.5*snv.width), par("usr")[1])),
       tick=F, line=0, labels="Nucleotides lost")
  axis(1, at=mean(c(snv.gap + (0.5*snv.width), par("usr")[2])),
       tick=F, line=0, labels="Nucleotides gained")

  # Add Y axis
  axis(2, at=log10(logscale.minor), labels=NA, tck=minor.tck, lwd=0.75)
  y.ax.at <- logscale.major[which(logscale.major >= 1)]
  clean.axis(2, at=log10(y.ax.at), labels=sapply(y.ax.at, clean.numeric.labels),
             title="Variant count", title.line=0.85)

  # Add legend
  text(x=0.5*snv.width, y=0.95*snv.k, labels="SNVs", font=3, pos=2,
       col=var.class.colors["snv"], xpd=T)
  text(x=-snv.gap-snv.width, y=max(del.xy$y),
       labels="Indels", font=3, pos=2, col=var.class.colors["indel"], xpd=T)
  text(x=-log10(50)-snv.gap-indel.gap-snv.width, y=max(loss.xy$y),
       labels="SVs", font=3, pos=2, col=var.class.colors["sv"], xpd=T)
  # n.balcpx <- sum(df[which(df$subclass %in% c("INV", "CPX", "CTX", "OTH")), -c(1:2)])
  # balcpx.lab <- clean.numeric.labels(n.balcpx)
  # text(x=par("usr")[2], y=0.8*par("usr")[4], col=annotation.color, cex=5/6, pos=2,
  #      labels=paste("Not shown:\n", balcpx.lab, "balanced &\ncomplex SVs"), xpd=T)
}

# Cowplot of SV sizes
sv.size.cowplot <- function(sv.sizes){
  ridgeplot(sv.sizes, xlims=log10(c(10, 5000000)), x.axis.side=NA,
            fill=sv.colors[names(sv.sizes)], bw.adj=sv.dens.bw[names(sv.sizes)],
            fancy.light.fill=adjust.brightness(sv.colors[names(sv.sizes)], d=1),
            fancy.median.color=adjust.brightness(sv.colors[names(sv.sizes)], d=-0.2),
            border=NA,
            parmar=c(2.2, 3.5, 0.1, 0.1))
  axis(1, at=log10(logscale.minor), tck=-0.01, labels=NA)
  clean.axis(1, at=log10(logscale.major.bp),
             labels=logscale.major.bp.labels[seq(1, length(logscale.major.bp), 2)],
             labels.at=log10(logscale.major.bp)[seq(1, length(logscale.major.bp), 2)],
             label.line=-0.9, title.line=0.2, title=bquote("SV size" ~ (log[10])))
}


# Plot AF distribution per variant class
plot.af.distribs <- function(af.df, breaks, colors=NULL, group.names=NULL, lwd=3,
                             y.title="Variant count", common.af=0.001,
                             parmar=c(2, 2.75, 0.25, 2.5)){
  # Prepare AF data
  af.dat <- lapply(1:nrow(af.df), function(i){
    step.function(x=breaks, y=unlist(log10(af.df[i, ])), offset=0)
  })

  # Get plot parameters
  xlims <- c(min(sapply(af.dat, function(l){l$x})) - 0.1,
             max(sapply(af.dat, function(l){l$x})))
  ylims <- c(max(c(0, (floor(10 * min(sapply(af.dat, function(l){l$y}))) / 10) - 0.1)),
             (ceiling(10 * max(sapply(af.dat, function(l){l$y}))) / 10) + 0.1)
  if(is.null(colors)){
    colors <- rev(greyscale.palette(length(af.dat)))
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar=parmar)

  # Add distinction between rare & common, if optioned
  if(!is.null(common.af)){
    abline(v=log10(common.af), col=annotation.color, lty=5)
    text(x=log10(common.af)+(0.035*diff(par("usr")[1:2])),
         y=par("usr")[4]-(0.04*diff(par("usr")[3:4])), cex=5/6, pos=2, labels="Rare")
    text(x=log10(common.af)-(0.035*diff(par("usr")[1:2])),
         y=par("usr")[4]-(0.04*diff(par("usr")[3:4])), cex=5/6, pos=4,
         labels="Common")
    text(x=log10(common.af)-(0.035*diff(par("usr")[1:2])),
         y=par("usr")[4]-(0.12*diff(par("usr")[3:4])), cex=5/6, pos=4,
         labels=paste("(AF>", round(100 * common.af, 1), "%)", sep=""))
  }

  # Add step functions and group names (if optioned)
  sapply(1:length(af.dat), function(i){
    points(af.dat[[i]]$x, af.dat[[i]]$y, col=colors[i], type="l", lwd=lwd)
  })

  # Add axes
  x.ax.at <- 0:floor(par("usr")[1])
  x.ax.labs <- head(c(1, 0.1, 0.01), length(x.ax.at))
  if(length(x.ax.at) > length(x.ax.labs)){
    x.ax.labs <- c(x.ax.labs, paste("10^", x.ax.at[-(1:length(x.ax.labs))], sep=""))
  }
  axis(1, at=log10(logscale.minor), tck=-0.01, lwd=0.75, labels=NA)
  clean.axis(1, at=x.ax.at, labels=x.ax.labs, parse.labels=TRUE, label.line=-0.9,
             title="Allele frequency (AF)", title.line=0)
  axis(2, at=log10(logscale.minor), tck=-0.01, lwd=0.75, labels=NA)
  clean.axis(2, at=log10(logscale.major.bp),
             labels=sapply(logscale.major.bp, clean.numeric.labels),
             title=y.title, title.line=0.85)

  # Add right Y-margin group labels
  if(!is.null(group.names)){
    yaxis.legend(group.names, x=par("usr")[2]-(0.05*diff(par("usr")[1:2])),
                 y.positions=sapply(af.dat, function(l){tail(l$y, 1)}),
                 sep.wex=0, min.label.spacing=0.1*diff(par("usr")[3:4]),
                 lower.limit=par("usr")[3]+0.025*(diff(par("usr")[3:4])),
                 upper.limit=par("usr")[4]-0.025*(diff(par("usr")[3:4])),
                 colors=NA, label.colors=colors, label.font=3)
  }
}

# Stacked barplot of AF bins by variant size
plot.size.by.af <- function(joint.d, bar.sep=0.1, parmar=c(2.6, 2.75, 0.25, 3.75)){
  # Convert data into proportions
  df <- joint.d$df
  af.bins <- joint.d$breaks
  size.bins <- sort(unique(clean.breaks(df$size)))
  names(size.bins) <- paste("ge", size.bins, sep="")
  plot.df <- do.call("cbind", lapply(size.bins, function(ge){
    k.by.af <- apply(df[which(df$size == paste("ge", ge, sep="")), -c(1:3)], 2, sum)
    k.by.af / sum(k.by.af)
  }))

  # Get plot parameters
  xlims <- c(0, length(size.bins))
  ylims <- c(0, 1)
  af.pal <- rev(viridis(length(af.bins), begin=0.05, end=0.95))

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar, xaxs="r")

  # Add bottom bin labels
  x.labs <- c("SNV", logscale.major.bp.labels[1:(length(size.bins)-1)])
  sapply(1:length(size.bins), function(x){
    text(x=x+0.1, y=-0.035, srt=40, pos=2, labels=x.labs[x], xpd=T, cex=5/6)
  })
  mtext(1, line=1.75, text="Variant size")

  # Add left Y axis
  clean.axis(2, at=seq(0, 1, 0.25), label.units="percent",
             title="Variant Proportion", title.line=0.6, label.line=-0.75)

  # Add right Y axis/legend hybrid
  af.labels.pct <- c("AF>10%", "AF<10%", "AF<1%", "AF<0.1%")
  if(length(af.bins) > length(af.labels.pct)){
    af.labels <- c(rep("", length(af.labels.pct)),
                   paste("AF<10^", -((length(af.labels.pct)+1):length(af.bins)), sep=""))
    parse.any <- TRUE
  }else{
    af.labels <- head(af.labels.pct, length(af.bins))
    parse.any <- FALSE
  }
  legend.y.pos <- (c(0, cumsum(plot.df[-nrow(plot.df), ncol(plot.df)])) + cumsum(plot.df[, ncol(plot.df)])) / 2
  legend.y.pos <- yaxis.legend(af.labels, x=ncol(plot.df)-bar.sep,
                               y.positions=legend.y.pos, parse.labels=parse.any,
                               upper.limit=0.975, colors=af.pal, sep.wex=0.5,
                               label.cex=5/6, min.label.spacing=0.1,
                               return.label.pos=TRUE)
  if(parse.any){
    text(x=ncol(plot.df)-bar.sep+0.5, y=rev(legend.y.pos)[1:length(af.labels.pct)],
         labels=af.labels.pct, cex=5/6, pos=4, xpd=T)
  }

  # Add bars last
  sapply(1:length(size.bins), function(x){
    rect(xleft=rep(x-1+bar.sep, nrow(plot.df)),
         xright=rep(x-bar.sep, nrow(plot.df)),
         ybottom=c(0, cumsum(plot.df[-nrow(plot.df), x])),
         ytop=cumsum(plot.df[, x]),
         col=af.pal, border=af.pal, lwd=0.5)
    rect(xleft=x-1+bar.sep, xright=x-bar.sep, ybottom=0, ytop=1, col=NA, xpd=T)
  })
}


###########
# RScript #
###########
# Parse command line arguments and opions
parser <- ArgumentParser(description="Plot variant site metric distributions")
parser$add_argument("--size-distrib", metavar=".tsv", type="character",
                    help="Precomputed binned variant size distribution")
parser$add_argument("--af-distrib", metavar=".tsv", type="character",
                    help="Precomputed binned variant AF distribution")
parser$add_argument("--joint-distrib", metavar=".tsv", type="character",
                    help="Precomputed variant size vs. AF distribution")
parser$add_argument("--sv-sites", metavar=".bed", type="character",
                    help="SV sites .bed file")
parser$add_argument("--common-af", metavar="float", default=0.01, type="numeric",
                    help="Allele frequency threshold for common variants")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    help="String or path to use as prefix for output plots",
                    default="./vcf_qc")
args <- parser$parse_args()

# # DEV:
# args <- list("size_distrib" = "~/scratch/gnomad_sites_test/gnomad.v4.1.size_distribution.merged.tsv.gz",
#              "af_distrib" = "~/scratch/gnomad_sites_test/gnomad.v4.1.af_distribution.merged.tsv.gz",
#              "joint_distrib" = "~/scratch/gnomad_sites_test/gnomad.v4.1.size_vs_af_distribution.merged.tsv.gz",
#              "sv_sites" = "~/scratch/YL.sv.site_metrics.dev.sv.sites.bed.gz",
#              "common_af" = 0.001,
#              "out_prefix" = "~/scratch/qc.test")

# Read distributions
size.d <- read.distrib(args$size_distrib)
af.d <- read.distrib(args$af_distrib)
joint.d <- read.distrib(args$joint_distrib, key.cols=1:3)
sv.sizes <- read.sv.sizes(args$sv_sites)

# Summary of variant counts
if(!is.null(size.d)){
  pdf(paste(args$out_prefix, "variant_count_bars.pdf", sep="."),
      height=2.25, width=2.85)
  plot.counts.by.vsc(size.d$df)
  dev.off()
}else if(!is.null(freq.d)){
  pdf(paste(args$out_prefix, "variant_count_bars.pdf", sep="."),
      height=2.25, width=2.85)
  plot.counts.by.vsc(freq.d$df)
  dev.off()
}

# Size plots
if(!is.null(size.d)){
  # Volcano of variant sizes
  pdf(paste(args$out_prefix, "variant_size_volcano.pdf", sep="."),
      height=2.25, width=4.1)
  plot.size.volcano(size.d, snv.gap=0.15, indel.gap=0)
  dev.off()
}
# Cowplot of SV sizes by subclass
if(!is.null(sv.sizes)){
  pdf(paste(args$out_prefix, "sv_size_cowplot.pdf", sep="."),
      height=2.25, width=2.5)
  sv.size.cowplot(sv.sizes)
  dev.off()
}



# AF plots
if(!is.null(af.d)){
  # Step function of AFs
  main.af.df <- do.call("rbind", lapply(names(var.class.colors), function(vc){
    apply(af.d$df[which(af.d$df$class == vc), -c(1:2)], 2, sum)
  }))
  pdf(paste(args$out_prefix, "af_distribs.pdf", sep="."),
      height=2.25, width=2.5)
  plot.af.distribs(main.af.df, breaks=log10(af.d$breaks), colors=var.class.colors,
                   group.names=paste(var.class.abbrev, "s", sep=""),
                   common.af=args$common_af)
  dev.off()

  # Supplementary step function of AFs for short variants by subclass
  short.af.df <- do.call("rbind", lapply(c("snv", "del", "ins"), function(vsc){
    apply(af.d$df[which(af.d$df$subclass == vsc), -c(1:2)], 2, sum)
  }))
  pdf(paste(args$out_prefix, "af_distribs.short_vars.pdf", sep="."),
      height=2.25, width=2.5)
  plot.af.distribs(short.af.df, breaks=log10(af.d$breaks),
                   colors=c(var.class.colors[["snv"]], indel.colors[c("del", "ins")]),
                   group.names=var.subclass.abbrevs[c("snv", "del", "ins")],
                   y.title="Short variant count", lwd=2, common.af=args$common_af)
  dev.off()

  # Supplementary step function of AFs for SVs by subclass
  svs.in.af <- intersect(names(sv.colors), af.d$df$subclass)
  sv.af.df <- do.call("rbind", lapply(svs.in.af, function(vsc){
    apply(af.d$df[which(af.d$df$subclass == vsc), -c(1:2)], 2, sum)
  }))
  pdf(paste(args$out_prefix, "af_distribs.svs.pdf", sep="."),
      height=2.25, width=2.5)
  plot.af.distribs(sv.af.df, breaks=log10(af.d$breaks),
                   colors=sv.colors[svs.in.af],
                   group.names=var.subclass.abbrevs[svs.in.af],
                   y.title="SV count", lwd=2, common.af=args$common_af)
  dev.off()
}

# 2D comparison of sizes vs. AFs
if(!is.null(joint.d)){
  pdf(paste(args$out_prefix, "size_vs_af_bars.pdf", sep="."),
      height=2.25, width=3.25)
  plot.size.by.af(joint.d)
  dev.off()
}

