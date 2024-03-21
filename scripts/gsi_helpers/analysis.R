library(dplyr)
library(readr)
library(stringr)
library(ggplot2)


# Read the TSV file with specified column names
df <- read.delim('/Users/noah/Desktop/DFCI_Data/gsi/data/gsc_dl.tsv', 
                 sep = '\t', 
                 col.names = c('patient', 'cancer', 'chrom', 'loc', 'germline_gene', 'somatic_gene','driver_likelihood'),
                 stringsAsFactors = FALSE)
gsc <- read.delim('/Users/noah/Desktop/DFCI_Data/gsi/data/VALab_germline_somatic_TSG_TSG.tsv', 
                  sep = '\t', 
                  stringsAsFactors = FALSE) %>% 
                  filter(germline_context == "coding" & somatic_context == "coding") %>% 
                  select('cancer','germline_gene','somatic_gene')

# Assuming gsc is your dataframe and cancer is the column
gsc$cancer <- str_to_title(gsc$cancer)
# Replace 'Renal' with 'Kidney' in gsc$cancer
gsc$cancer <- gsub("Renal", "Kidney", gsc$cancer)


calculate_fisher_odds <- function(df, cancer_type, germline_gene, somatic_gene) {
  # Filter dataframe by cancer type
  if (!is.na(cancer_type) && cancer_type != "All") {
    df_filtered <- df %>% filter(cancer == cancer_type)
  } else {
    df_filtered <- df
  }
  
  df_gene <- df_filtered

  # Function to create contingency table for a patient
  create_cont_table <- function(patient_data) {
    germ <- ifelse(germline_gene %in% patient_data$germline_gene, 1, 0)
    som <- ifelse(somatic_gene %in% patient_data$somatic_gene, 1, 0)
    return(data.frame(germline_gene = germ, somatic_gene = som))
  }
  
  # Group data by patient and create contingency tables for each patient
  patient_cont_tables <- lapply(split(df_gene, df_gene$patient), create_cont_table)
  
  # Combine all patient contingency tables into one
  combined_cont_table <- do.call(rbind, patient_cont_tables)
  
  # Create final contingency table
  final_cont_table <- table(combined_cont_table$germ, combined_cont_table$som)
  
  # Label rows and columns
  colnames(final_cont_table) <- c("Somatic_0", "Somatic_1")
  rownames(final_cont_table) <- c("Germline_0", "Germline_1")
  
  # Print contingency table if needed
  #print("Contingency Table:")
  print(cancer_type)
  print(final_cont_table)
  
  # Perform Fisher's exact test
  fisher_result <- fisher.test(final_cont_table)
  
  # Calculate Odds Ratio and 95% Confidence Interval
  odds_ratio <- fisher_result$estimate
  conf_interval <- fisher_result$conf.int
  
  # Print Fisher's exact test result if needed
  # print("Fisher's Exact Test:")
  # print(fisher_result)
  
  # Return p-value, Odds Ratio, and 95% Confidence Interval
  return(list(
    fisher_result$p.value,
    fisher_result$estimate,
    fisher_result$conf.int
  ))
}



p_value <- c()
OR <- c()
counter <- 0
file_conn <- file("/Users/noah/Desktop/output.txt", open = "w")  # "a" for append, "w" to overwrite
close(file_conn)
file_conn <- file("/Users/noah/Desktop/output.txt", open = "a")  # "a" for append

for (i in 1:nrow(gsc)) {
  row <- gsc[i, ]
  cancer_type <- row$cancer
  germline_gene_name <- row$germline_gene
  somatic_gene_name <- row$somatic_gene
  #We have to make sure that both the Germline and Somatic genes are there at 
  #some point (although they don't need to be at the same site/patient) 
  
  # Filter the dataframe for the current cancer type
  df_filtered <- df %>% filter(cancer == cancer_type)
  
  # Check if germline_gene_name is present in germline_gene column
  germline_present <- germline_gene_name %in% df_filtered$germline_gene
  
  # Check if somatic_gene_name is present in somatic_gene column
  somatic_present <- somatic_gene_name %in% df_filtered$somatic_gene
  
  if (germline_present & somatic_present) {
    result <- calculate_fisher_odds(df, cancer_type = cancer_type, germline_gene = germline_gene_name, somatic_gene = somatic_gene_name)
    p_value <- c(p_value,result[[1]])
    OR <- c(OR,result[[2]])
    #writeLines(paste(cancer_type,germline_gene_name,somatic_gene_name,log2(result[[2]]),-log10(result[[1]])),file_conn)
    #print(result)
  }
}
close(file_conn)


##################Graphing#####################
create_volcano_plot <- function(p_values, odds_ratios) {
  # Convert p-values to -log10(p-values)
  neg_log_p <- -log10(p_values)
  
  # Calculate log2(Odds Ratio)
  log2_odds_ratio <- log2(odds_ratios)
  
  # Create dataframe for plotting
  volcano_df <- data.frame(
    log2_OR = log2_odds_ratio,
    neg_log_p = neg_log_p
  )
  
  # Create volcano plot
  volcano_plot <- ggplot(volcano_df, aes(x = log2_OR, y = neg_log_p)) +
    geom_point(size = 2) +
    labs(
      x = "log2(Odds Ratio)",
      y = "-log10(p-value)",
      title = "Volcano Plot Coding-Coding TSG-TSG Convergence"
    ) +
    theme_minimal() +
    scale_x_continuous(limits = c(-5, 15)) +  # Center x-axis at 0 with range -4 to 4
    scale_y_continuous(limits = c(0, max(neg_log_p) + 1))  # Start y-axis at 0
  
  # Print plot
  print(volcano_plot)
}

# Example data
#p_values <- c(1.5e-16, 3e-9, 1.75e-68)
#odds_ratios <- c(Inf, 41, 96)

# Call the function to create volcano plot
create_volcano_plot(p_value, OR)
