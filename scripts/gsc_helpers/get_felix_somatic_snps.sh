#!/bin/bash

# Download the zip file
#wget https://www.science.org/doi/suppl/10.1126/science.abg5601/suppl_file/science.abg5601_tables_s2_to_s23.zip

# Unzip the file
#unzip science.abg5601_tables_s2_to_s23.zip
#rm combined_data.tsv combined_data.tmp

#get somatic noncoding portions
gsc_table="/Users/noah/Desktop/DFCI_Data/gsi/data/VALab_germline_somatic.tsv"
cut -f1,4,5 "$gsc_table" | grep 'noncoding' | cut -f2 > somatic_noncoding.list

# Define the Python script to extract sheet names and process each sheet
read -r -d '' PYTHON_SCRIPT <<'EOF'
import pandas as pd

# Load the Excel file
excel_file = pd.ExcelFile('science.abg5601_tables_s2_to_s23/science.abg5601_tables_s2_to_s23.xlsx')

# Loop through each sheet in the Excel file
for sheet_name in excel_file.sheet_names:
    cancer_type = sheet_name.split()[-1]
    cancer_type_set = set(['Breast','Colorectal','Prostate','Lung','Kidney'])
    if cancer_type not in cancer_type_set:
        continue
    # Read the data from the sheet into a DataFrame
    df = pd.read_excel(excel_file, sheet_name=sheet_name)
    # Export the data to a TSV file
    df.to_csv(f'{sheet_name.split()[-1]}.tsv', sep='\t', index=False)
    #print(sheet_name.split()[-1])
EOF

# Execute the Python script
python3 -c "$PYTHON_SCRIPT"

# Combine all TSV files into one
cat *.tsv | grep -E 'Breast|Colorectal|Prostate|Kidney|Lung' \
| grep -v 'Table' \
| awk -F'\t' '$3 != "coding"' \
| cut -f1,2,7 | awk -F'\t' -v OFS='\t' '{ split($3, chr, "-"); print chr[1],chr[2], $1, $2 }' \
| awk -F'\t' -v OFS='\t' '{ split($1, a, ":"); print a[1],a[2], $2, $3, $4 }' \
| grep -Fwf somatic_noncoding.list | sort -n > combined_data.tmp

#rm somatic_noncoding.list

# Remove individual TSV files
rm Breast.tsv Prostate.tsv Lung.tsv Kidney.tsv Colorectal.tsv

mv combined_data.tmp gsc_noncoding_snps_hg38.bed

