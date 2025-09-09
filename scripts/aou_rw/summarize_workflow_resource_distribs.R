#!/usr/bin/env Rscript

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Summarize distribution of compute resource usage per Cromwell task
# Note: this is a helper script called within summarize_workflow_resources.sh

options(scipen=1000, stringsAsFactors=F)

# Load task-level resources
res <- read.table(commandArgs(trailingOnly=TRUE)[1],
                  header=T, sep="\t", comment.char="", check.names=F)
colnames(res)[1] <- gsub("#", "", colnames(res)[1])

# Define strata to summarize
strata <- unique(res[, c("task", "resource")])

# Compute summary metrics per stratum
out.df <- as.data.frame(do.call("rbind", lapply(1:nrow(strata), function(i){
  sub <- res[which(res$task == strata$task[i]
            & res$resource == strata$resource[i]), ]
  n.tasks <- nrow(sub)
  if(all(is.na(sub$allocated))){
    alloc <- NA
  }else{
    alloc <- mean(sub$allocated, na.rm=T)
  }
  c(sub$task[1], sub$resource[1], n.tasks, alloc, as.numeric(summary(sub$peak_used)))
})))
colnames(out.df) <- c("task", "resource", "n_calls", "mean_allocated", "min_used",
                      "q1_used", "median_used", "mean_used", "q3_used", "max_used")
out.df <- out.df[with(out.df, order(task, resource)), ]

# Report to stdout
write.table(out.df, file="", sep="\t", row.names=F, quote=F)
