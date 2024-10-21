#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Plot cohort-wide intake QC metrics


#########
# Setup #
#########
# Load necessary libraries and constants
options(scipen=1000, stringsAsFactors=F)
require(argparse, quietly=TRUE)
require(G2C, quietly=TRUE)
require(survival, quietly=TRUE)
require(vioplot, quietly=TRUE)
G2C::load.constants("all")


##################
# Data Functions #
##################
# Load & clean QC dataframe
load.qc.df <- function(tsv.in){
  qc.df <- read.table(tsv.in, check.names=F, header=T, sep="\t", comment.char="")
  colnames(qc.df)[1] <- gsub("^#", "", colnames(qc.df)[1])
  rownames(qc.df) <- qc.df$G2C_id
  qc.df$G2C_id <- NULL
  qc.df$cohort_type <- remap(qc.df$cohort, cohort.type.map)
  qc.df$single_cancer <- qc.df$cancer
  qc.df$single_cancer[grepl(";", qc.df$single_cancer, fixed=T)] <- "multiple"
  bool.cols <- grep("_qc_pass", colnames(qc.df))
  if(length(bool.cols) > 0){
    qc.df[, bool.cols] <- apply(qc.df[, bool.cols], 2, remap,
                                map=c("True" = TRUE, "False" = FALSE))
  }
  return(qc.df)
}


######################
# Plotting Functions #
######################
# Grafpop margin bar
grafpop.margin.bar <- function(qc.df, ordered.pops=NULL, label.cex=5.5/6){
  prep.plot.area(c(0, 1), c(0, nrow(qc.df)), parmar=c(2.5, 0.25, 0.5, 7))
  if(is.null(ordered.pops)){
    ordered.pops <- intersect(names(pop.colors), unique(qc.df$intake_qc_pop))
  }
  pop.k.df <- do.call("rbind", lapply(ordered.pops, function(pop){
    data.frame("pop"=pop, rev(sort(table(qc.df$intake_qc_subpop[which(qc.df$intake_qc_pop == pop)]))))
  }))
  colnames(pop.k.df) <- c("pop", "subpop", "n")
  pop.k.df$cumsum <- cumsum(pop.k.df$n)
  super.k.df <- data.frame("n" = sapply(ordered.pops, function(pop){
    length(which(qc.df$intake_qc_pop == pop))
  }))
  super.k.df$pop <- rownames(super.k.df)
  super.k.df$cumsum <- cumsum(super.k.df$n)
  rect(xleft=0, xright=1, ybottom=c(0, pop.k.df$cumsum[-nrow(pop.k.df)]),
       ytop=pop.k.df$cumsum, col=subpop.colors[as.character(pop.k.df$subpop)],
       border=NA, bty="n")
  rect(xleft=0.75, xright=1, ybottom=c(0, super.k.df$cumsum[-nrow(super.k.df)]),
       ytop=super.k.df$cumsum, col=pop.colors[ordered.pops], border="white", lwd=2)
  pop.pcts <- round(100 * super.k.df$n / sum(super.k.df$n), 1)
  legend.labels <- paste(pop.names.short[ordered.pops], " (", pop.pcts, "%)", sep="")
  yaxis.legend(legend.labels, x=1,
               y.positions=(c(0, super.k.df$cumsum[-nrow(super.k.df)]) + super.k.df$cumsum)/2,
               sep.wex=0.25, min.label.spacing=0.09*diff(par("usr")[3:4]),
               label.cex=label.cex,
               colors=pop.colors[ordered.pops])
}

# Sex ploidy margin bar
sex.margin.bar <- function(qc.df, label.cex=5.5/6){
  sex.df <- data.frame(row.names = c("female", "other", "male"),
                       "sex" = c("female", "other", "male"),
                       "n" = c(length(which(qc.df$inferred_sex == "female")),
                               length(which(qc.df$inferred_sex == "other")),
                               length(which(qc.df$inferred_sex == "male"))))
  sex.df$cumsum <- cumsum(sex.df$n)
  sex.legend.labels <- c(
    paste("XX (", round(100 * sex.df["female", "n"] / nrow(qc.df), 1), "%)", sep=""),
    paste("Other (", round(100 * sex.df["other", "n"] / nrow(qc.df), 1), "%)", sep=""),
    paste("XY (", round(100 * sex.df["male", "n"] / nrow(qc.df), 1), "%)", sep="")
  )
  prep.plot.area(c(0, 1), c(0, nrow(qc.df)), parmar=c(2.5, 0.25, 0.5, 7))
  rect(xleft=0, xright=1, ybottom=c(0, sex.df$cumsum[-nrow(sex.df)]),
       ytop=sex.df$cumsum, col=sex.colors[rownames(sex.df)], border=NA, bty="n")
  yaxis.legend(sex.legend.labels, x=1,
               y.positions=(c(0, sex.df$cumsum[-nrow(sex.df)]) + sex.df$cumsum)/2,
               sep.wex=0.25, min.label.spacing=0.09*diff(par("usr")[3:4]),
               label.cex=label.cex, colors=sex.colors[rownames(sex.df)])
}

# Swarms of all autosomal ploidy estimates
plot.autosomal.ploidy <- function(qc.df, parmar=c(1.9, 2.1, 1, 0.5)){
  # Ensure libraries are loaded
  require(vioplot)

  # Get plot parameters
  ploidy.df <- qc.df[, grep("^chr[0-9]+_ploidy$", colnames(qc.df))]
  min.ploidy <- min(c(0.75, min(ploidy.df, na.rm=T)))
  max.ploidy <- max(c(3.25, max(ploidy.df, na.rm=T)))

  # Prepare plot area
  prep.plot.area(c(0, 22), c(min.ploidy, max.ploidy), parmar=parmar)
  abline(h=c(1, 2, 3), lty=c(5, 1, 5), col="gray80")
  aneu.y.adj <- 0.05 * diff(par("usr")[3:4])
  text(x=par("usr")[1], y=c(1, 3) + (aneu.y.adj * c(1, -1)), pos=4,
       labels=c("Monosomy", "Trisomy"), cex=5/6, font=3, col="gray80")
  mtext(1, text="Chromosome", line=1)
  clean.axis(2, title="Ploidy", infinite=TRUE, title.line=0.2)
  mtext(3, text="Autosomal ploidy estimates")

  # Add points for each contig
  sapply(1:22, function(x){
    # Add X axis label
    axis(1, at=x-0.5, tick=F, line=-0.9, cex.axis=5/6, labels=x)

    # Get all ploidy values
    pl.vals <- as.numeric(ploidy.df[, paste("chr", x, "_ploidy", sep="")])

    # For sake of plotting, don't plot a point for every single sample
    # Instead, approximate the "normal" range (Q1/3 ± 3IQR) as a vioplot
    # Then add beeswarm for outlier points
    pl.out.bool <- label.outliers(pl.vals)
    pl.in <- pl.vals[which(!pl.out.bool)]
    pl.out <- pl.vals[which(pl.out.bool)]
    vioplot(pl.in, add=T, at=x-0.5, drawRect=F, pchMed=NA, bty="n",
            col=DFCI.colors[["paleblue"]])
    pl.med <- median(pl.vals, na.rm=T)
    segments(x0=x-0.5-0.125, x1=x-0.5+0.125, y0=pl.med, y1=pl.med, pch=3)
    points(x=rep(x-0.5, length(which(pl.out.bool))), y=pl.out, cex=0.125, xpd=T)
  })
}

# Plot a matrix/grid of selected features vs. one categorical variable (cohort|cancer)
# TODO: need to implement this for batch as the primary column
qc.grid <- function(qc.df, primary.variable){
  # Format primary and secondary variables
  if(primary.variable == "cancer"){
    prim <- cancer <- as.character(explode.by.cancer(qc.df$cancer)$cancer)
    cohort <- as.character(explode.by.cancer(qc.df$cancer, qc.df$simple_cohort)$x)
    sex <- as.character(explode.by.cancer(qc.df$cancer, qc.df$inferred_sex)$x)
    pop <- as.character(explode.by.cancer(qc.df$cancer, qc.df$intake_qc_pop)$x)
    age <- as.numeric(explode.by.cancer(qc.df$cancer, qc.df$age)$x)
    coverage <- as.numeric(explode.by.cancer(qc.df$cancer, qc.df$mean_coverage)$x)
    wgd <- as.numeric(explode.by.cancer(qc.df$cancer, qc.df$wgd_score)$x)
  }else{
    cancer <- as.character(qc.df$single_cancer)
    prim <- cohort <- as.character(qc.df$simple_cohort)
    sex <- as.character(qc.df$inferred_sex)
    pop <- as.character(qc.df$intake_qc_pop)
    age <- as.numeric(qc.df$age)
    coverage <- as.numeric(qc.df$mean_coverage)
    wgd <- as.numeric(qc.df$wgd_score)
  }
  prim.order <- names(sort(table(prim)))
  prim.titles <- list("cancer" = cancer.names, "cohort" = cohort.names.short)[[primary.variable]]
  sec.order <- setdiff(c("cancer", "cohort", "sex", "pop", "age", "coverage"),
                       primary.variable)

  # Configure oscillating greyscale for cohort colors on the fly
  cohort.k <- sort(table(qc.df$simple_cohort), decreasing=TRUE)
  cohort.colors <- greyscale.palette(length(cohort.k), oscillate=TRUE)
  names(cohort.colors) <- names(cohort.k)

  # Prep plot layout
  n.rows <- length(prim.order)
  n.cols <- length(sec.order) + 1
  layout(matrix(1:n.cols, nrow=1), widths=c(2, rep(1.4, n.rows-1)))
  top.y.margin <- 2.5
  inner.margin <- 0.3

  # First panel is always simple barplot of total abundance per row + right margin labels
  if(primary.variable == "cancer"){
    prim.colors <- cancer.colors[prim.order]
    names(prim.colors) <- prim.titles[prim.order]
  }else if(primary.variable == "cohort"){
    prim.colors <- cohort.colors
    names(prim.colors) <- prim.titles[names(cohort.colors)]
  }
  stacked.barplot(prim.titles[prim], color=prim.colors,
                  outer.borders=prim.colors[prim.titles[prim.order]],
                  x.axis.side=NA, orient="left",
                  custom.major.order=prim.titles[prim.order],
                  add.legend=FALSE, annotate.counts=TRUE,
                  parmar=c(0.1, 3, top.y.margin, 4.5))
  axis(3, at=mean(par("usr")[1:2]), tick=F, line=-0.8, labels="Genomes")

  for(sec in sec.order){
    # Get values for each secondary variable
    if(sec == "cancer"){
      s.names <- cancer.names
      s.v <- s.names[cancer]
      s.col <- cancer.colors
      names(s.col) <- s.names[names(s.col)]
      s.title <- "Cancer type"
      plot.type <- "bar"
      sort.minor <- TRUE

    }else if(sec == "cohort"){
      s.names <- cohort.names.short
      s.v <- s.names[cohort]
      s.k <- sort(table(s.v), decreasing=TRUE)
      s.col <- cohort.colors
      names(s.col) <- s.names[names(cohort.colors)]
      s.title <- "Cohort"
      plot.type <- "bar"
      sort.minor <- TRUE

    }else if(sec == "sex"){
      s.names <- genetic.sex.names
      s.v <- s.names[sex]
      s.col <- sex.colors
      names(s.col) <- s.names[names(s.col)]
      s.title <- "Genetic sex"
      plot.type <- "bar"
      sort.minor <- FALSE

    }else if(sec == "pop"){
      s.names <- names(pop.colors)
      names(s.names) <- s.names
      s.v <- s.names[pop]
      s.col <- pop.colors
      names(s.col) <- s.names[names(s.col)]
      s.title <- "Genetic\nancestry"
      plot.type <- "bar"
      sort.minor <- TRUE

    }else if(sec == "age"){
      s.v <- lapply(prim.order, function(p){
        p.a <- as.numeric(age[which(prim == p)])
        p.a <- p.a[which(!is.na(p.a))]
      })
      xlims <- NULL
      s.title <- "Age (years)"
      plot.type <- "ridge"

    }else if(sec == "coverage"){
      s.v <- lapply(prim.order, function(p){
        p.c <- as.numeric(coverage[which(prim == p)])
        p.c <- p.c[which(!is.na(p.c))]
      })
      s.title <- "WGS coverage (x)"
      xlims <- c(0, min(c(80, quantile(unlist(s.v), probs=0.999))))
      plot.type <- "ridge"

    }

    if(plot.type == "bar"){
      stacked.barplot(major.values=prim, minor.values=s.v, colors=s.col,
                      as.proportion=T, add.legend=F, add.major.label=F,
                      minor.labels.on.bars=TRUE, x.axis.side=NA,
                      minor.label.letter.width=0.07, sort.minor=sort.minor,
                      parmar=c(0.1, inner.margin, top.y.margin, inner.margin))

      if(sec == "sex"){
        abline(v=0.5, lty=3)
        clean.axis(3, at=0.5, labels="50:50", title=s.title, title.line=0.3,
                   label.line=-0.65)
      }else{
        axis(3, at=0.5, tick=F, line=-0.5, labels=s.title)
      }
    }else if(plot.type == "ridge"){
      if(primary.variable == "cancer"){
        ridge.fill <- cancer.colors[prim.order]
        ridge.light <- sapply(cancer.palettes[prim.order], function(v){v["light1"]})
        ridge.border <- sapply(cancer.palettes[prim.order], function(v){v["dark2"]})
        median.color <- sapply(cancer.palettes[prim.order], function(v){v["light2"]})
      }else if(primary.variable == "cohort"){
        ridge.fill <- cohort.colors[prim.order]
        ridge.light <- adjust.brightness(ridge.fill, 0.15)
        ridge.border <- adjust.brightness(ridge.fill, -0.3)
        median.color <- adjust.brightness(ridge.fill, 0.3)
      }
      ridgeplot(s.v, x.title=s.title, x.axis.side=3, xlims=xlims, y.axis=FALSE,
                bw.adj=0.5, yaxs="i", hill.overlap=-0.1, hill.bottom=0.1,
                border.lwd=1, fill=ridge.fill, fancy.light.fill=ridge.light,
                border=ridge.border, fancy.median.color=median.color,
                parmar=c(0.1, inner.margin+0.1, top.y.margin, inner.margin+0.1))
    }else{
      prep.plot.area(0:1, 0:1, rep(0, 4))
    }
  }
}


###########
# RScript #
###########
# Parse command line arguments and options
parser <- ArgumentParser(description="Visualize intake QC metrics")
parser$add_argument("--qc-tsv", metavar=".tsv", type="character", required=TRUE,
                    help="Intake QC .tsv")
parser$add_argument("--pass-column", metavar="string", type="character",
                    help=paste("Intake QC .tsv will be filtered according to",
                               "boolean (or boolean-coercible) 'true' in",
                               "this column prior to plotting. Column values",
                               "are not case sensitive. Can be specified",
                               "multiple times"), action="append")
parser$add_argument("--out-prefix", metavar="path", type="character",
                    default="G2C.intake_qc",
                    help="String or path to use as prefix for output plots")
args <- parser$parse_args()

# # DEV:
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.intake_qc.non_aou.post_qc_batching.tsv.gz",
#              "pass_column" = c("global_qc_pass", "batch_qc_pass"),
#              "out_prefix" = "~/scratch/dfci-g2c.intake_qc.local_test")

# Load data
qc.df <- load.qc.df(args$qc_tsv)


# Set shared barplot parameters
bar.parmars <- list(c(0, 4.35, 1.25, 2.75), c(0, 4.35, 2, 2.75))
fail.bar.pdf.dims <- c(4.5, 4.5)
cancer.key.cols <- cancer.colors
names(cancer.key.cols) <- cancer.names[names(cancer.colors)]
cohort.k <- sort(table(qc.df$simple_cohort), decreasing=TRUE)
cohort.cols <- greyscale.palette(length(cohort.k), oscillate=TRUE)
names(cohort.cols) <- cohort.names.short[names(cohort.k)]
cohort.key.cols <- cohort.cols


# Handle --pass-column behavior, if optioned
if(!is.null(args$pass_column) & length(args$pass_column) > 0){
  # Collect list of row indexes failing any of --pass-column
  fail.rows <- c()
  for(c.name in args$pass_column){
    if(c.name %in% colnames(qc.df)){
      if(is.numeric(qc.df[, c.name])){
        fail.rows <- c(fail.rows, which(!as.logical(qc.df[, c.name])))
      }else{
        fail.rows <- c(fail.rows, which(!as.logical(toupper(as.character(qc.df[, c.name])))))
      }
    }else{
      stop(paste("Column name '", c.name, "', specified as --pass-column,",
                 "but this column could not be located in --qc-tsv. Exiting."))
    }
  }
  fail.rows <- sort(unique(fail.rows))
  fail.df <- qc.df[fail.rows, ]

  # Before excluding samples, first generate barplots of failure rate per cohort & cancer
  # Barplot by cancer type, colored by cohort
  pdf(paste(args$out_prefix, "failure_rate_per_cancer", "pdf", sep="."),
      height=fail.bar.pdf.dims[1], width=fail.bar.pdf.dims[2])
    stacked.barplot(cancer.names[explode.by.cancer(fail.df$cancer)$cancer],
                    cohort.names.short[explode.by.cancer(fail.df$cancer, fail.df$simple_cohort)$x],
                    colors=cohort.cols, x.title="Samples failing QC",
                    x.label.line=-0.85, x.axis.tck = -0.0125,
                    x.title.line=0, annotate.counts=TRUE, add.legend=FALSE,
                    major.legend=TRUE, major.legend.colors=cancer.key.cols, major.legend.xadj=-0.02,
                    minor.labels.on.bars=TRUE, minor.label.letter.width=0.05,
                    minor.label.cex=4/6, parmar=bar.parmars[[2]])
  dev.off()
  # Barplot by cohort, colored by cancer
  pdf(paste(args$out_prefix, "failure_rate_per_cohort", "pdf", sep="."),
      height=fail.bar.pdf.dims[1], width=fail.bar.pdf.dims[2])
  stacked.barplot(cohort.names.short[explode.by.cancer(fail.df$cancer, fail.df$simple_cohort)$x],
                  cancer.names[explode.by.cancer(fail.df$cancer)$cancer],
                  colors=cancer.key.cols, x.title="Samples failing QC",
                  x.label.line=-0.85, x.axis.tck = -0.0125,
                  x.title.line=0, annotate.counts=TRUE, add.legend=FALSE,
                  major.legend=TRUE, major.legend.colors=cohort.cols, major.legend.xadj=-0.02,
                  minor.labels.on.bars=TRUE, minor.label.letter.width=0.05, sort.minor=TRUE,
                  minor.label.cex=4/6, parmar=bar.parmars[[2]])
  dev.off()



  # Lastly, drop samples outright from qc.df before continuing
  qc.df <- qc.df[-fail.rows, ]
}


# Reset some barplot parameters after excluding QC failures
cancer.k <- sort(table(qc.df$single_cancer), decreasing=TRUE)
cancer.key.cols <- cancer.key.cols[cancer.names[names(cancer.k)]]
cohort.k <- sort(table(qc.df$simple_cohort), decreasing=TRUE)
cohort.cols <- greyscale.palette(length(cohort.k), oscillate=TRUE)
names(cohort.cols) <- cohort.names.short[names(cohort.k)]
cohort.type.k <- sort(table(qc.df$cohort_type), decreasing=TRUE)
cohort.type.cols <- greyscale.palette(length(cohort.type.k), oscillate=TRUE)
names(cohort.type.cols) <- cohort.type.names.short[names(cohort.type.k)]


# Grafpop coordinates by ancestry
apply(t(combn(1:3, 2)), 1, function(gd.idxs){
  idx.1 <- gd.idxs[1]
  idx.2 <- gd.idxs[2]
  gp.out.prefix <- paste(args$out_prefix, ".grafpop_gd", idx.1, "_gd", idx.2,
                      "_by_subpop", sep="")
  for(device in c("pdf", "png")){
    if(device == "pdf"){
      pdf(paste(gp.out.prefix, "pdf", sep="."), height=2.5, width=4)
    }else if(device == "png"){
      png(paste(gp.out.prefix, "png", sep="."),
          height=2.5*300, width=4*300, res=300, family="Arial")
    }
    layout(matrix(1:2, nrow=1), widths=c(7, 5))
    scatterplot(qc.df[, paste("grafpop_GD", idx.1, sep="")],
                qc.df[, paste("grafpop_GD", idx.2, sep="")],
                title="Genetic ancestry",
                subpop.colors[qc.df$intake_qc_subpop],
                x.title=paste("Genetic distance", idx.1),
                y.title=paste("Genetic distance", idx.2),
                parmar=c(2.25, 2.5, 1, 0.35),
                x.label.line=-0.7, y.label.line=-0.6, x.title.line=0.25)
    pop.order <- names(sort(sapply(unique(qc.df$intake_qc_pop), function(pop){
      mean(qc.df[which(qc.df$intake_qc_pop == pop),
                 paste("grafpop_GD", idx.2, sep="")], na.rm=T)
    })))
    grafpop.margin.bar(qc.df, pop.order)
    dev.off()
  }
})


# Visualize sex ploidy distributions
for(device in c("pdf", "png")){
  if(device == "pdf"){
    pdf(paste(args$out_prefix, "sex_ploidy.pdf", sep="."), height=2.5, width=4)
  }else if(device == "png"){
    png(paste(args$out_prefix, "sex_ploidy.png", sep="."),
        height=2.5*300, width=4*300, res=300, family="Arial")
  }
  layout(matrix(1:2, nrow=1), widths=c(7, 5))
  scatterplot(qc.df$chrX_ploidy, qc.df$chrY_ploidy,
              colors=sex.colors[qc.df$inferred_sex],
              title="Genetic sex",
              x.title="chrX ploidy", y.title="chrY ploidy",
              x.title.line=0.15,
              parmar=c(2.1, 2.6, 1, 0.5),
              x.label.line=-0.7, y.label.line=-0.6)
  sex.margin.bar(qc.df)
  dev.off()
}


# Autosomal ploidy estimates
for(device in c("pdf", "png")){
  if(device == "pdf"){
    pdf(paste(args$out_prefix, "autosomal_ploidy.pdf", sep="."),
        height=2.25, width=5)
  }else if(device == "png"){
    png(paste(args$out_prefix, "autosomal_ploidy.png", sep="."),
        height=2.25*300, width=5*300, res=300, family="Arial")
  }
  plot.autosomal.ploidy(qc.df)
  dev.off()
}


# Numeric metrics with no special coloring
hist.pdf.dims <- c(1.7, 2.3)
for(metric.info in list(c("hq_hom", "High-qual. hom. GTs", "count"),
                        c("hq_hom_rate", "High-qual. hom. GT rate", "percent"),
                        c("hq_het", "High-qual. het. GTs", "count"),
                        c("hq_het_rate", "High-qual. het. GT rate", "percent"),
                        c("charr", "Contamination", "percent"),
                        c("mean_ref_ab_hom_alt", "Hom. GT allele balance", "percent"),
                        c("heterozygosity_rate", "Heterozygosity", "percent"),
                        c("inconsistent_ab_het_rate", "Bad het. allele balance", "percent"),
                        c("mean_coverage", "Mean coverage", "other"),
                        c("median_coverage", "Median coverage", "other"),
                        c("wgd_score", "Dosage bias", "other"),
                        c("pct_genome_nondiploid", "Nondiploid genome frac.", "percent"),
                        c("age", "Age at intake", "count"),
                        c("years_to_last_contact", "Years of follow-up", "count"))){
  vals <- as.numeric(qc.df[, metric.info[1]])
  vals <- vals[which(!is.na(vals) & !is.infinite(vals))]
  title <- as.character(metric.info[2])
  pdf(paste(args$out_prefix, metric.info[1], "pdf", sep="."),
      height=hist.pdf.dims[1], width=hist.pdf.dims[2])
  density.w.outliers(vals, title=title, x.label.units=metric.info[3],
                     bw.adj=2, style="hist", col=DFCI.colors[["paleblue"]],
                     outlier.lwd=0.5, parmar=c(1.1, 2.5, 1, 1))
  dev.off()
}


# Format data for main barplots
bar.pdf.dims <- c(3.5, 3)
cancer.panel.height.ratio <- c(1/5, 1)
cohort.panel.height.ratio <- c(1/4, 1)
non.cancer.phenos <- c("control", "unknown")
cancer.bar.subdfs <- list(qc.df[which(qc.df$cancer %in% non.cancer.phenos), ],
                          qc.df[which(!qc.df$cancer %in% non.cancer.phenos), ])
largest.cohort <- names(cohort.k)[1]
cohort.bar.subdfs <- list(qc.df[which(qc.df$simple_cohort == largest.cohort), ],
                          qc.df[which(qc.df$simple_cohort != largest.cohort), ])
cohort.type.k <- sort(table(qc.df$cohort_type), decreasing=TRUE)
cohort.type.cols <- greyscale.palette(length(cohort.type.k), oscillate=TRUE)
names(cohort.type.cols) <- cohort.type.names.short[names(cohort.type.k)]
all.cohorts <- names(sort(table(qc.df$cohort), decreasing=TRUE))
cohort.key.cols <- cohort.type.cols[cohort.type.names.short[cohort.type.map[all.cohorts]]]
names(cohort.key.cols) <- cohort.names.short[all.cohorts]
cohort.key.cols["Other"] <- cohort.type.cols["Oncology"]
frac.multi <- sum(grepl(";", qc.df$cancer)) / sum(qc.df$cancer %in% c("control", "unknown"))
annot.color <- "gray70"


# Barplot of samples per cancer, colored by cohort
pdf(paste(args$out_prefix, "cohort_contributions_per_cancer", "pdf", sep="."),
    height=bar.pdf.dims[1], width=bar.pdf.dims[2])
layout(matrix(2:1, nrow=2, ncol=1), heights=rev(cancer.panel.height.ratio))
sapply(1:2, function(s){
  subdf <- cancer.bar.subdfs[[s]]
  stacked.barplot(cancer.names[explode.by.cancer(subdf$cancer)$cancer],
                  cohort.type.names.short[explode.by.cancer(subdf$cancer, subdf$cohort_type)$x],
                  colors=cohort.type.cols,
                  x.title=if(s==2){"Genomes per cancer type"}else{""},
                  x.axis.tck=-0.025/cancer.panel.height.ratio[s], x.label.line=-0.85,
                  x.title.line=0, annotate.counts=TRUE, add.legend=FALSE,
                  major.legend=TRUE, major.legend.colors=cancer.key.cols,
                  minor.labels.on.bars=TRUE, minor.label.letter.width=0.05,
                  minor.label.cex=4/6, parmar=bar.parmars[[s]])
  if(s==2){
    axis(4, at=par("usr")[3]-1.25, tick=F, las=2, hadj=1, line=1, cex.axis=5/6,
         labels=paste(round(100 * frac.multi, 1), "% of cases have\nmultiple cancers"),
         xpd=T, col.axis=annot.color)
    axis(4, at=par("usr")[3]-4.25, tick=F, las=2, hadj=1, line=1, cex.axis=5/6,
         labels=paste("Mean = ", prettyNum(round(mean(cancer.k), 0), big.mark=","),
                      "\nper cancer type", sep=""), xpd=T, col.axis=annot.color)
  }
})
dev.off()


# Barplot of samples by cohort, colored by cancer
pdf(paste(args$out_prefix, "cancers_per_cohort", "pdf", sep="."),
    height=bar.pdf.dims[1], width=bar.pdf.dims[2])
layout(matrix(1:2, nrow=2, ncol=1), heights=cohort.panel.height.ratio)
sapply(1:2, function(s){
  subdf <- cohort.bar.subdfs[[s]]
  stacked.barplot(cohort.names.short[subdf$simple_cohort],
                  cancer.names[subdf$cancer], colors=cancer.key.cols, inner.borders=cancer.key.cols,
                  x.title=if(s==1){"Genomes per cohort"}else{""},
                  x.axis.tck=-0.025/cohort.panel.height.ratio[s], x.label.line=-0.85,
                  x.title.line=0, annotate.counts=TRUE, add.legend=FALSE,
                  major.legend=TRUE, major.legend.colors=cohort.key.cols,
                  minor.labels.on.bars=TRUE, minor.label.letter.width=0.07,
                  minor.label.cex=4/6, custom.minor.order=names(cancer.key.cols),
                  parmar=bar.parmars[[if(s==1){2}else{1}]])
  if(s==2){
    axis(4, at=par("usr")[3]-0.75, tick=F, las=2, hadj=1, line=1, cex.axis=5/6,
         labels=paste("N = ", prettyNum(nrow(qc.df), big.mark=","), "\n",
                      "total genomes", sep=""), xpd=T, col.axis=annot.color)
  }
})
dev.off()


# Barplot of stage for cancer cases
# Note: uses hist.pdf.dims (this is intentional)
pdf(paste(args$out_prefix, "cancer_stages", "pdf", sep="."),
    height=hist.pdf.dims[1], width=hist.pdf.dims[2]-0.3)
stacked.barplot(stage.names[qc.df$stage], colors=stage.colors,
                x.title="Stage at Dx",
                add.legend=F, custom.major.order=stage.names,
                x.label.line=-0.8, x.title.line=0.1,
                parmar=c(0.5, 3, 2.5, 0.75))
dev.off()


# Kaplan-Meier curves for all cases versus all controls
surv.models <- lapply(c("IV", "III", "II", "I"), function(stage){
  summary(survfit(Surv(years_left_censored, years_to_last_contact, abs(vital_status-1)) ~ 1,
                  data=qc.df[which(qc.df$stage == stage), ]))
})
names(surv.models) <- c("IV", "III", "II", "I")
surv.models[["no_stage"]] <- summary(survfit(Surv(years_left_censored, years_to_last_contact, abs(vital_status-1)) ~ 1,
                                            data=qc.df[which(qc.df$cancer != "control" & !is.na(qc.df$stage)), ]))
surv.models[["control"]] <- summary(survfit(Surv(years_left_censored, years_to_last_contact, abs(vital_status-1)) ~ 1,
                              data=qc.df[which(qc.df$cancer == "control"), ]))
# Note: uses hist.pdf.dims (this is intentional)
pdf(paste(args$out_prefix, "km_surv_by_stage", "pdf", sep="."),
    height=hist.pdf.dims[1], width=hist.pdf.dims[2]+0.3)
km.curve(surv.models,
         colors=c(stage.colors[c("IV", "III", "II", "I")], cancer.colors[c("pancan", "control")]),
         ci.alpha=0, title="Overall survival", y.title="", time.is.days=FALSE,
         legend.label.spacing=0.15, legend.label.cex=5/6,
         legend.names=c(stage.names.long[names(surv.models)[1:4]], "Stage N.R.", "Controls"),
         xlim=c(0, 8), x.tck=-0.025, x.label.line=-1, x.title.line=-0.25,
         km.lwd=4, parmar=c(1.75, 1.25, 1, 5))
dev.off()


# Summary grid of key metrics by cancer type & cohort
sapply(c("cancer", "cohort"), function(primary.variable){
  pdf(paste(args$out_prefix, "qc_grid", primary.variable, "pdf", sep="."),
      height=if(primary.variable == "cancer"){4}else{2.5}, width=7.2)
  qc.grid(qc.df, primary.variable)
  dev.off()
})


# TODO: implement qc matrix for all batches (full page layout for supplement)
# Also split in ~half for double-wide layout on powerpoint slide
