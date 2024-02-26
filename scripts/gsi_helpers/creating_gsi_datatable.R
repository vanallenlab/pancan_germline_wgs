library(dplyr)
library(readr)

match_df <- read_tsv("/Users/noah/Desktop/HMF.combined_metadata_and_manifest.tsv") %>% select(hmfPatientId,sampleId)
all_samples_df <- read_tsv("/Users/noah/Desktop/HMF.hg19_germline_vcf_processing.terra_manifest.tsv")

# Extract filename from the URI in the 'germline_vcf' column
all_samples_df$filename <- basename(all_samples_df$germline_vcf)
# Create new column 'germlineVCFs' with the desired path
all_samples_df$germlineVCFs <- paste0('gs://hmf-dr-355-us-central1/hg19_germline_vcfs/', all_samples_df$filename)
all_samples_df$germlineVCF_index <- paste0('gs://hmf-dr-355-us-central1/hg19_germline_vcfs/', all_samples_df$filename,'.tbi')

# Optionally, you can remove the 'filename' column if it's no longer needed
all_samples_df <- all_samples_df[, -which(names(all_samples_df) == "filename")]
merged_df <- left_join(all_samples_df,match_df)

merged_df$somaticDrivers <- paste0('gs://fc-9cb68074-23c1-4bb3-9ef2-7363efd1fb40/somatic_tarball/',merged_df$sampleId,'/purple/',merged_df$sampleId,'.purple.driver.catalog.somatic.tsv')


df_processed <- merged_df %>%
  group_by(hmfPatientId,primaryTumorLocation,primaryTumorSubLocation,primaryTumorType,
           g2c_cancer_type,germlineVCFs,germlineVCF_index) %>%
  summarise(    somaticDrivers = sprintf("[%s]", paste0('"', paste(unique(somaticDrivers), collapse = '","'), '"')),
                sampleId = sprintf("[%s]", paste0('"', paste(unique(sampleId), collapse = '","'), '"')))

#merged_df <- merged_df %>% select(sampleId,everything())
write_tsv(df_processed,"/Users/noah/Desktop/hmf_gsi_cohort.tsv")
