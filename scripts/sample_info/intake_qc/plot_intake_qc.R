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
# args <- list("qc_tsv" = "~/scratch/dfci-g2c.intake_qc.merged.test.tsv",
#              "pass_column" = NULL,
#              "out_prefix" = "~/scratch/dfci-g2c.intake_qc.local_test")

# Load data
qc.df <- load.qc.df(args$qc_tsv)

# Filter data to --pass-column, if optioned
if(!is.null(args$pass_column) & length(args$pass_column) > 0){
  for(c.name in args$pass_column){
    if(c.name %in% colnames(qc.df)){
      if(is.numeric(qc.df[, c.name])){
        qc.df <- qc.df[which(as.logical(qc.df[, c.name])), ]
      }else{
        qc.df <- qc.df[which(as.logical(toupper(as.character(qc.df[, c.name])))), ]
      }
    }else{
      stop(paste("Column name '", c.name, "', specified as --pass-column,",
                 "but this column could not be located in --qc-tsv. Exiting."))
    }
  }
}

# Grafpop coordinates by ancestry
apply(t(combn(1:3, 2)), 1, function(gd.idxs){
  idx.1 <- gd.idxs[1]
  idx.2 <- gd.idxs[2]
  pdf(paste(args$out_prefix, ".grafpop_gd", idx.1, "_gd", idx.2,
            "_by_subpop.pdf", sep=""),
      height=2.5, width=4)
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
})


# Visualize sex ploidy distributions
pdf(paste(args$out_prefix, "sex_ploidy.pdf", sep="."),
    height=2.5, width=4)
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


# Autosomal ploidy estimates
pdf(paste(args$out_prefix, "autosomal_ploidy.pdf", sep="."),
    height=2.25, width=5)
plot.autosomal.ploidy(qc.df)
dev.off()


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


# Set shared barplot parameters
bar.parmars <- list(c(0.25, 4.25, 1.25, 2.75), c(0, 4.25, 2, 2.75))
bar.pdf.dims <- c(3.5, 3)
cancer.panel.height.ratio <- c(1/4.65, 1)
cohort.panel.height.ratio <- c(1/4, 1)
cancer.key.cols <- cancer.colors
non.cancer.phenos <- c("control", "unknown")
cancer.bar.subdfs <- list(qc.df[which(qc.df$cancer %in% non.cancer.phenos), ],
                          qc.df[which(!qc.df$cancer %in% non.cancer.phenos), ])
names(cancer.key.cols) <- cancer.names[names(cancer.colors)]
cohort.k <- sort(table(qc.df$simple_cohort), decreasing=TRUE)
largest.cohort <- names(cohort.k)[1]
cohort.bar.subdfs <- list(qc.df[which(qc.df$simple_cohort == largest.cohort), ],
                          qc.df[which(qc.df$simple_cohort != largest.cohort), ])
cohort.key.cols <- greyscale.palette(length(cohort.k), oscillate=TRUE)
names(cohort.key.cols) <- cohort.names.short[names(cohort.k)]

# Barplot of samples per cancer, colored by cohort
pdf(paste(args$out_prefix, "cohort_contributions_per_cancer", "pdf", sep="."),
    height=bar.pdf.dims[1], width=bar.pdf.dims[2])
layout(matrix(2:1, nrow=2, ncol=1), heights=rev(cancer.panel.height.ratio))
sapply(1:2, function(s){
  subdf <- cancer.bar.subdfs[[s]]
  stacked.barplot(cancer.names[explode.by.cancer(subdf$cancer)$cancer],
                  cohort.names.short[explode.by.cancer(subdf$cancer, subdf$simple_cohort)$x],
                  colors=cohort.key.cols,
                  x.title=if(s==2){"Genomes per cancer type"}else{""},
                  x.axis.tck=-0.025/cancer.panel.height.ratio[s], x.label.line=-0.85,
                  x.title.line=0, annotate.counts=TRUE, add.legend=FALSE,
                  major.legend=TRUE, major.legend.colors=cancer.key.cols,
                  minor.labels.on.bars=TRUE, minor.label.letter.width=0.06,
                  minor.label.cex=4/6, parmar=bar.parmars[[s]])
  if(s==2){
    axis(4, at=par("usr")[3]-1.25, tick=F, las=2, hadj=1, line=1, cex.axis=5/6,
         labels=paste("N = ", prettyNum(nrow(qc.df), big.mark=","), "\n",
                      "total genomes", sep=""), xpd=T)
  }
})
dev.off()


# Barplot of samples by cohort, colored by cancer
pdf(paste(args$out_prefix, "cancers_per_cohort", "pdf", sep="."),
    height=bar.pdf.dims[1], width=bar.pdf.dims[2])
layout(matrix(1:2, nrow=2, ncol=1), heights=cohort.panel.height.ratio)
sapply(1:2, function(s){
  subdf <- cohort.bar.subdfs[[s]]
  subdf$cancer[grepl(";", subdf$cancer, fixed=T)] <- "multiple"
  stacked.barplot(cohort.names.short[subdf$simple_cohort],
                  cancer.names[subdf$cancer], colors=cancer.key.cols,
                  x.title=if(s==1){"Genomes per cohort"}else{""},
                  x.axis.tck=-0.025/cohort.panel.height.ratio[s], x.label.line=-0.85,
                  x.title.line=0, annotate.counts=TRUE, add.legend=FALSE,
                  major.legend=TRUE, major.legend.colors=cohort.key.cols,
                  minor.labels.on.bars=TRUE, minor.label.letter.width=0.07,
                  minor.label.cex=4/6, sort.minor=TRUE,
                  parmar=bar.parmars[[if(s==1){2}else{1}]])
  # if(s==2){
  #   axis(4, at=par("usr")[3]-0.75, tick=F, las=2, hadj=1, line=1, cex.axis=5/6,
  #        labels=paste("N = ", prettyNum(nrow(qc.df), big.mark=","), "\n",
  #                     "total genomes", sep=""), xpd=T)
  # }
})
dev.off()


# Barplot of stage for cancer cases
# Note: uses hist.pdf.dims (this is intentional)
pdf(paste(args$out_prefix, "cancer_stages", "pdf", sep="."),
    height=hist.pdf.dims[1], width=hist.pdf.dims[2]-0.3)
stacked.barplot(stage.names[qc.df$stage], colors=stage.colors,
                x.title="Stage at Dx",
                add.legend=F, custom.order=stage.names,
                x.label.line=-0.8, x.title.line=0.1)
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

