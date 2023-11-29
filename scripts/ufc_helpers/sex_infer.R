library(dplyr)

estimate_chromosomes_processing <- function(file_path) {
  # Create an R data frame with chromosome numbers and lengths
  chromosome_data <- data.frame(
    Chromosome = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10",
                   "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
                   "chr21", "chr22", "chrX", "chrY"),
    Length = c(248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636,
               138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345,
               81282511, 78077248, 59128983, 63025520, 48129895, 50735781, 156040895, 57227415)
  )
  
  # Read the TSV file into a data frame
  df <- read.table(file_path, header = FALSE, sep = "\t", col.names = c("Chromosome", "V2", "V3", "V4", "V5", "V6"), fill = TRUE)
  
  # Remove the 5th and 6th columns
  df <- df[, -c(5, 6)]
  
  # Find the row index where the specified header pattern is located
  header_row <- grep("CONTIG", df$Chromosome)
  
  # Omit rows before the header row (if header is found)
  df <- df[(header_row + 1):nrow(df), ]
  
  # Group by V1 and calculate the sum of V4 for each group
  df <- df %>%
    group_by(Chromosome) %>%
    summarize(Sum_V4 = sum(as.numeric(V4)))
  
  # Merge the two data frames based on the "Chromosome" column
  df <- merge(df, chromosome_data, by = "Chromosome")
  
  # Calculate normalized read coverage
  df$norm_read_cov <- df$Sum_V4 / df$Length
  
  # Return the merged data frame with normalized read coverage
  return(df)
}

estimate_chromosomes <- function(df) {
  # Extract 'norm_read_cov' values for autosomes and sex chromosomes
  autosomal_cov <- df$norm_read_cov[!df$Chromosome %in% c("chrX", "chrY")]
  chrX_cov <- df$norm_read_cov[df$Chromosome == "chrX"]
  chrY_cov <- df$norm_read_cov[df$Chromosome == "chrY"]
  
  # Calculate mean and standard deviation of autosomal counts
  autosomal_mean <- mean(autosomal_cov)
  autosomal_sd <- sd(autosomal_cov)
  
  # Set thresholds based on mean and standard deviation
  threshold_one_copy <- autosomal_mean - 3 * autosomal_sd
  threshold_two_copies <- autosomal_mean - 2 * autosomal_sd
  threshold_three_copies <- autosomal_mean + 2 * autosomal_sd
  
  # Initialize estimates
  estimated_chrX <- 0
  estimated_chrY <- 0
  
  # Iterate through each threshold and add 1 copy for every threshold it passes
  if (max(chrX_cov) > threshold_one_copy) {
    estimated_chrX <- estimated_chrX + 1
  }
  
  if (max(chrX_cov) > threshold_two_copies) {
    estimated_chrX <- estimated_chrX + 1
  }
  
  if (max(chrX_cov) > threshold_three_copies) {
    estimated_chrX <- estimated_chrX + 1
  }
  
  if (max(chrY_cov) > threshold_one_copy) {
    estimated_chrY <- estimated_chrY + 1
  }
  
  if (max(chrY_cov) > threshold_two_copies) {
    estimated_chrY <- estimated_chrY + 1
  }
  
  if (max(chrY_cov) > threshold_three_copies) {
    estimated_chrY <- estimated_chrY + 1
  }
  
  # Return the estimates in the specified format
  return(paste(estimated_chrX, estimated_chrY, sep = "-"))
}


# Extract the file path from the command-line arguments
file_path <- commandArgs(trailingOnly = TRUE)[1]

# Example usage
result_df <- estimate_chromosomes_processing(file_path)
result <- estimate_chromosomes(result_df)
cat(result)
