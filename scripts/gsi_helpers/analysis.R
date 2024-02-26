library(dplyr)
library(readr)
library(stringr)

# Reading dataframe from a CSV file (replace 'your_dataframe_file.csv' with the actual file path)
# Specify the column names
column_names <- c('patient', 'cancer', 'chrom', 'loc', 'germ', 'som')

# Read the TSV file with specified column names
df <- read.delim('/Users/noah/Desktop/pancan_germline_wgs/scripts/gsi_helpers/data/gsc_filtered.tsv', 
                 sep = '\t', 
                 col.names = column_names,
                 stringsAsFactors = FALSE)
df$som <- str_trim(df$som, side = "left")

df_germ <- df %>% filter(germ == "APC" | germ == "VHL" | germ == "BRCA1")

df_som <- df %>% filter(som == "APC" | som == "VHL" | som == "BRCA1")
