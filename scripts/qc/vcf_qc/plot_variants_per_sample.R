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
require(DescTools, quietly=TRUE)
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

# Filter GT counts based on column values
filter.gt.counts <- function(gt.counts, vc=NULL, vsc=NULL, freq=NULL){
  # Subset by class
  if(!is.null(vc)){
    gt.counts <- gt.counts[which(gt.counts$class %in% vc), ]
  }

  # Subset by subclass
  if(!is.null(vsc)){
    if(is.na(vsc)){
      gt.counts <- gt.counts[which(is.na(gt.counts$subclass)), ]
    }else{
      gt.counts <- gt.counts[which(gt.counts$subclass %in% vsc), ]
    }
  }

  # Subset by frequency
  if(!is.null(freq)){
    gt.counts <- gt.counts[which(gt.counts$freq_bin %in% freq), ]
  }

  gt.counts
}

# Compute a summary metric across a subset of gt counts
calc.gt.ss <- function(gt.counts, samples, zyg=NA, summary.fx=median){
  # Sum counts for each sample
  k <- sapply(samples, function(sid){
    s.k <- gt.counts[which(gt.counts$sample == sid),
                     if(is.na(zyg)){c("het", "hom")}else{zyg}]
    if(length(s.k) == 0){0}else{sum(as.numeric(unlist(s.k)))}
  })
  if(is.null(summary.fx)){
    return(k)
  }else{
    summary.fx(k)
  }
}

# Main wrapper function to collect all per-sample summary stats
gather.count.sumstats <- function(gt.counts, samples, pop=NULL, pheno=NULL){
  as.data.frame(do.call("rbind", lapply(unique(gt.counts$class), function(vc){
    vc.prefix <- paste("variants_per_genome", vc, sep=".")
    vc.df <- filter.gt.counts(gt.counts, vc=vc)

    do.call("rbind", lapply(unique(vc.df$subclass), function(vsc){
      vsc.prefix <- paste(vc.prefix, vsc, sep=".")
      vsc.df <- filter.gt.counts(vc.df, vsc=vsc)

      do.call("rbind", lapply(unique(vsc.df$freq_bin), function(freq){
        freq.prefix <- paste(vsc.prefix, paste("AF", freq, sep="_"), sep=".")
        freq.df <- filter.gt.counts(vsc.df, freq=freq)

        do.call("rbind", lapply(c(NA, "het", "hom"), function(zyg){
          zyg.prefix <- sub(".NA", "", paste(freq.prefix, zyg, sep="."), fixed=T)
          med.v <- calc.gt.ss(freq.df, samples=samples, zyg=zyg, summary.fx=median)
          phi.v <- calc.gt.ss(freq.df, samples=samples, zyg=zyg,
                              summary.fx=RLCtools::dynamic.range)
          zyg.sub.res <- data.frame("analysis"=zyg.prefix,
                                    "measure"=c("median", "dynamic_range"),
                                    "value"=c(med.v, phi.v), "n"=length(samples))

          if(!is.null(pop)){
            zyg.sub.add <- do.call("rbind", lapply(unique(pop), function(p){
              pop.prefix <- sub(".NA", "", paste(zyg.prefix, p, sep="."), fixed=T)
              pop.samples <- names(pop)[which(pop == p)]
              med.v <- calc.gt.ss(freq.df, samples=pop.samples, zyg=zyg, summary.fx=median)
              phi.v <- calc.gt.ss(freq.df, samples=pop.samples, zyg=zyg,
                                  summary.fx=RLCtools::dynamic.range)
              zyg.sub.add <- data.frame("analysis"=pop.prefix,
                                        "measure"=c("median", "dynamic_range"),
                                        "value"=c(med.v, phi.v),
                                        "n"=length(pop.samples))

              if(!is.null(pheno)){
                zyg.sub.add.pw <- do.call("rbind", lapply(unique(pheno), function(ph){
                  pw.prefix <- sub(".NA", "", paste(pop.prefix, ph, sep="."), fixed=T)
                  pheno.samples <- names(pheno)[which(pheno == ph)]
                  pw.samples <- intersect(pop.samples, pheno.samples)
                  if(length(pw.samples) == 0){
                    return(NULL)
                  }
                  med.v <- calc.gt.ss(freq.df, samples=pw.samples, zyg=zyg, summary.fx=median)
                  phi.v <- calc.gt.ss(freq.df, samples=pw.samples, zyg=zyg,
                                      summary.fx=RLCtools::dynamic.range)
                  data.frame("analysis"=pw.prefix, "measure"=c("median", "dynamic_range"),
                             "value"=c(med.v, phi.v), "n"=length(pw.samples))
                }))
                zyg.sub.add <- rbind(zyg.sub.add, zyg.sub.add.pw)
              }
              return(zyg.sub.add)
            }))
            zyg.sub.res <- rbind(zyg.sub.res, zyg.sub.add)
          }

          if(!is.null(pheno)){
            zyg.sub.add <- do.call("rbind", lapply(unique(pheno), function(p){
              pheno.prefix <- sub(".NA", "", paste(zyg.prefix, p, sep="."), fixed=T)
              pheno.samples <- names(pheno)[which(pheno == p)]
              med.v <- calc.gt.ss(freq.df, samples=pheno.samples, zyg=zyg, summary.fx=median)
              phi.v <- calc.gt.ss(freq.df, samples=pheno.samples, zyg=zyg,
                                  summary.fx=RLCtools::dynamic.range)
              data.frame("analysis"=pheno.prefix, "measure"=c("median", "dynamic_range"),
                         "value"=c(med.v, phi.v), "n"=length(pheno.samples))
            }))
            zyg.sub.res <- rbind(zyg.sub.res, zyg.sub.add)
          }
          return(zyg.sub.res)
        }))
      }))
    }))
  })))
}


######################
# Plotting Functions #
######################
# Helper function to plot a single count waterfall; called within main.waterfall()
plot.count.waterfall <- function(gt.counts, vc, pop=NULL, pheno=NULL,
                                 samples.sorted=NULL, pop.space.wex=0.025,
                                 pheno.space.wex=0.005,
                                 parmar=c(0.25, 3.25, 0.5, 0.15)){
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
  total.per.sample <- apply(counts, 1, sum, na.rm=T)

  # Get plot parameters
  pop.spacer <- n.samples * pop.space.wex
  pheno.spacer <- n.samples * pheno.space.wex
  xlims <- c(-0.01*n.samples,
             n.samples + (pop.spacer * length(pop.breaks.x)) + (pheno.spacer * length(pheno.breaks.x)))
  ylims <- c(0, max(total.per.sample))
  hom.col <- adjust.color.hsb(var.class.colors[vc], s=0.01, b=-0.05)
  het.col <- adjust.color.hsb(var.class.colors[vc], s=-0.01, b=0.05)

  # Determine x position for each sample
  prev.x <- 0
  sample.xleft <- c()
  for(i in 1:n.samples){
    if(i %in% pop.breaks.x){
      prev.x <- prev.x + pop.spacer
    }else if(i %in% pheno.breaks.x){
      prev.x <- prev.x + pheno.spacer
    }
    sample.xleft[samples.sorted[i]] <- prev.x
    prev.x <- prev.x + 1
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar)

  # Add bars
  sapply(samples.sorted, function(sid){
    rect(xleft=sample.xleft[sid],
         xright=sample.xleft[sid]+1,
         ybottom=c(0, counts[sid, 1]),
         ytop=cumsum(as.numeric(counts[sid, 1:2])),
         col=c(hom.col, het.col), border=c(hom.col, het.col), lwd=0.25)
  })

  # Add groupwise medians
  unique.s.feats <- unique(s.feats[, c("pop", "pheno")])
  sapply(1:nrow(unique.s.feats), function(i){
    pop <- unique.s.feats$pop[i]
    pheno <- unique.s.feats$pheno[i]
    sids <- rownames(s.feats)[which(s.feats$pop == pop & s.feats$pheno == pheno)]
    group.df <- data.frame("x" = sample.xleft[sids]+0.5,
                           "k" = total.per.sample[sids],
                           row.names=sids)
    group.df <- group.df[order(group.df$x), ]
    lab.x.at <- median(group.df$x)
    lab.y.at <- max(group.df$k[which(group.df$x >= quantile(group.df$x, 0.25)
                      & group.df$x <= quantile(group.df$x, 0.75))])
    group.k.med <- median(group.df$k)
    text(x=lab.x.at, y=lab.y.at, cex=5/6, xpd=T,
         col=MixColor(var.class.colors[vc], "black", 2/3),
         labels=clean.numeric.labels(round(group.k.med, 0),
                                     acceptable.decimals=if(group.k.med <= 2000){0}else{1}))
  })

  # Add axis
  clean.axis(2, title=paste(var.class.abbrevs[vc], "s", sep=""),
             label.units="count", infinite.positive=T, max.ticks=4,
             title.line=0.75, cex.title=7.5/6, cex.axis=1)
}

# Helper function to add groupwise markers to waterfall plot
add.waterfall.markers <- function(order.df, pop.map, pheno.map,
                                  pop.space.wex=0.025, pheno.space.wex=0.005,
                                  group.labels=c("Ancestry", "Pheno."),
                                  rect.hex.buffer=0.125,
                                  parmar=c(1.35, 3.25, 0.25, 0.15)){
  # Determine breaks between groups
  s.feats <- data.frame("order" = 1:nrow(order.df), row.names=rownames(order.df))
  s.feats$pop <- if(is.null(pop)){NA}else{pop[rownames(order.df)]}
  s.feats$pheno <- if(is.null(pop)){NA}else{pheno[rownames(order.df)]}
  n.samples <- nrow(s.feats)
  pop.spacer <- pop.space.wex * n.samples
  pheno.spacer <- pheno.space.wex * n.samples
  last.x <- 0
  pop.labs <- pheno.labs <- pop.starts <- pop.ends <- pheno.starts <- pheno.ends <- c()
  for(r in 2:n.samples){
    last.x <- last.x + 1
    if(s.feats$pop[r] != s.feats$pop[r-1]){
      pop.ends <- c(pop.ends, last.x)
      pheno.ends <- c(pheno.ends, last.x)
      last.x <- last.x + pop.spacer
      pop.starts <- c(pop.starts, last.x)
      pheno.starts <- c(pheno.starts, last.x)
      pop.labs <- c(pop.labs, s.feats$pop[r-1])
      pheno.labs <- c(pheno.labs, s.feats$pheno[r-1])
    }else{
      if(s.feats$pheno[r] != s.feats$pheno[r-1]){
        pheno.ends <- c(pheno.ends, last.x)
        last.x <- last.x + pheno.spacer
        pheno.starts <- c(pheno.starts, last.x)
        pheno.labs <- c(pheno.labs, s.feats$pheno[r-1])
      }
    }
    if(r == n.samples){
      if(length(pop.ends) > 0){
        pop.ends <- c(pop.ends, last.x)
        pop.labs <- c(pop.labs, s.feats$pop[r])
      }
      if(length(pheno.ends) > 0){
        pheno.ends <- c(pheno.ends, last.x)
        pheno.labs <- c(pheno.labs, s.feats$pheno[r])
      }
    }
  }
  do.pop <- length(pop.starts) > 0
  if(do.pop){
    pop.starts <- c(0, pop.starts)
  }
  do.pheno <- length(pheno.starts) > 0
  if(do.pheno){
    pheno.starts <- c(0, pheno.starts)
  }

  # Get plot parameters
  xlims <- c(-0.01*n.samples, max(c(pop.ends, pheno.ends, n.samples)))
  ylims <- c(0, ncol(order.df)-1)
  if(all(pop.labs %in% names(pop.colors))){
    pop.pal <- pop.colors[names(pop.map)]
    pop.names <- pop.names.short[names(pop.map)]
  }else{
    pop.pal <- categorical.rainbow(length(pop.map))
    pop.names <- names(pop.map)
    names(pop.names) <- names(pop.pal) <- names(pop.map)
  }
  if(all(pheno.labs %in% names(cancer.colors))){
    pheno.pal <- cancer.colors[names(pheno.map)]
    pheno.names <- cancer.names[names(pheno.map)]
  }else{
    pheno.pal <- greyscale.palette(length(pheno.map)+2)[2:(length(pheno.map)+1)]
    pheno.names <- names(pheno.map)
    names(pheno.names) <- names(pheno.pal) <- names(pheno.map)
  }

  # Prep plot area
  prep.plot.area(xlims, ylims, parmar)

  # Draw rectangles with labels
  y.start <- 0
  if(do.pop){
    rect(xleft=pop.starts, xright=pop.ends,
         ybottom=y.start + rect.hex.buffer,
         ytop=y.start + 1 - rect.hex.buffer,
         col=pop.pal[pop.labs], border=pop.pal[pop.labs])
    pop.rect.df <- data.frame("start" = pop.starts, "end" = pop.ends,
                              "len" = pop.ends-pop.starts, "lab" = pop.labs)
    pop.rect.df <- pop.rect.df[order(-pop.rect.df$len), ]
    sapply(unique(pop.labs), function(pop){
      pop.lab.x <- mean(as.numeric(pop.rect.df[which(pop.rect.df$lab == pop),
                                               c("start", "end")]))
      text(x=pop.lab.x, y=y.start+0.5, labels=pop.names[pop],
           col=optimize.label.color(pop.pal[pop], cutoff=0.8))
    })
    y.start <- y.start + 1
  }
  if(do.pheno){
    rect(xleft=pheno.starts, xright=pheno.ends,
         ybottom=y.start + rect.hex.buffer, ytop=y.start + 1 - rect.hex.buffer,
         col=pheno.pal[pheno.labs], border=pheno.pal[pheno.labs])
    pheno.rect.df <- data.frame("start" = pheno.starts, "end" = pheno.ends,
                                "len" = pheno.ends-pheno.starts, "lab" = pheno.labs)
    pheno.rect.df <- pheno.rect.df[order(-pheno.rect.df$len), ]
    sapply(unique(pheno.labs), function(pheno){
      pheno.lab.x <- mean(as.numeric(head(pheno.rect.df[which(pheno.rect.df$lab == pheno),
                                                        c("start", "end")], 1)))
      text(x=pheno.lab.x, y=y.start+0.5, labels=pheno.names[pheno],
           col=optimize.label.color(pheno.pal[pheno], cutoff=0.8))
    })
  }

  # Add labels
  sapply(1:length(group.labels), function(y){
    axis(2, at=y-0.5, tick=F, line=-1.2, labels=group.labels[y], las=2)
  })
}

# Wrapper for main waterfall plot
main.waterfall <- function(gt.counts, out.prefix, pop=NULL, pheno=NULL,
                           pop.space.wex=0.025, pheno.space.wex=0.005){
  # Determine number of panels and figure sizing
  vcs <- intersect(names(var.class.names), unique(gt.counts$class))
  n.panels <- length(vcs)
  has.pop <- !is.null(pop)
  has.pheno <- !is.null(pheno)
  bottom.tracks <- sum(c(has.pop, has.pheno))
  has.bottom.tracks <- bottom.tracks > 0
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
  order.df <- order.df[do.call(order, order.df), ]
  samples.sorted <- rownames(order.df)

  # Define plot parameters
  # Target dimensions for 3 rows + 2 marker fields: 6.6" wide x 4" tall
  top.to.bottom.ratio <- 2.5
  height.scalar <- 4 / ((3*top.to.bottom.ratio) + 1)
  pdf.width <- 6.6
  pdf.height <- (top.to.bottom.ratio*height.scalar*n.panels) + (bottom.tracks*height.scalar/2)

  # Generate waterfall plot
  needs.legend <- TRUE
  pdf(paste(out.prefix, "variants_per_genome.waterfall.pdf", sep="."),
      height=pdf.height, width=pdf.width)
  layout(matrix(1:(n.panels+has.bottom.tracks), byrow=T, ncol=1),
         heights=c(rep(top.to.bottom.ratio, n.panels), if(has.bottom.tracks){1}))
  for(vc in vcs){
    waterfall.parmar <- c(0.25, 3.25, 0.5, 0.15)
    if(!has.bottom.tracks & vc == vcs[n.panels]){
      parmar[1] <- 1.5
    }
    plot.count.waterfall(gt.counts, vc, pop, pheno, samples.sorted,
                         pop.space.wex, pheno.space.wex, parmar=waterfall.parmar)
    if(needs.legend){
      legend("topright", fill=c("gray80", "gray60"), cex=7/6,
             legend=c("Heterozygous", "Homozygous"), bty="n", border=NA)
      needs.legend <- FALSE
    }
  }
  if(has.bottom.tracks){
    add.waterfall.markers(order.df, pop.rank, pheno.rank,
                          pop.space.wex, pheno.space.wex)
  }
  mtext(paste(prettyNum(length(samples), big.mark=","), "germline genomes"),
        1, line=0.25, cex=5/6)
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

# Collect median counts per sample by class, subclass, frequency, zygosity, pop, and pheno
count.ss <- gather.count.sumstats(gt.counts, samples, pop, pheno)

# TODO: maybe add Sankey diagrams here

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
if(nrow(ic.ss) > 0){
  count.ss <- as.data.frame(rbind(count.ss, ic.ss))
}

# Write summary statistics to output file
colnames(count.ss)[1] <- paste("#", colnames(count.ss)[1], sep="")
write.table(count.ss, paste(args$out_prefix, "summary_metrics.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
