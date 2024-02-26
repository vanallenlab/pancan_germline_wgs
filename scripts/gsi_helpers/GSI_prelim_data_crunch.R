library(dplyr)
library(readr)
library(ggplot2)

df1 <- read_tsv("/Users/noah/Desktop/pancan_germline_wgs/scripts/gsi_helpers/data/hiv_all_samples.tsv")
df2 <- read_tsv("/Users/noah/Desktop/pancan_germline_wgs/scripts/gsi_helpers/data/warnings_and_failures.tsv")
gene_list <- read_tsv("/Users/noah/Desktop/pancan_germline_wgs/scripts/gsi_helpers/data/germline_genes.tsv")
df1$gene_loc <- paste(df1$CHROM,":",df1$POSITION,sep="")


df3 <- merge(df1,df2,by="gene_loc")
df4 <- df3 %>% filter(flag == "good" | flag == "flag")
df5 <- df4 %>% filter(GENE %in% gene_list$germline_gene) %>% select(Patient,Cancer,GENE,CHROM,POSITION)
#write_tsv(df5,"/Users/noah/Desktop/pancan_germline_wgs/scripts/gsi_helpers/data/filtered_germline_patient_data.tsv")
# Group by Patient and count occurrences
df_counts <- df5 %>%
  group_by(Patient) %>%
  summarise(count = n())

# Create a histogram of the counts
ggplot(df_counts, aes(x = count)) +
  geom_histogram(binwidth = 1, fill = "blue", color = "black") +
  labs(title = "Histogram of Patient Counts", x = "Count", y = "Frequency")
