import pandas as pd

def create_bed_file(input_file, cancer_type, output_file):
    # Read the input CSV file into a DataFrame
    df = pd.read_csv(input_file, sep='\t')
    
    # Rename columns to standardized BED format
    df = df.rename(columns={"CHR_ID": "chrom", "CHR_POS": "pos", "MAPPED_GENE": "gene"})
    
    # Prepare the BED columns
    bed_df = pd.DataFrame({
        "chrom": "chr" + df["chrom"].astype(str),
        "start": df["pos"] - 1,  # BED is 0-based
        "end": df["pos"],
        "gene": df["gene"].str.split(',').str[0],
        "risk_allele": df["STRONGEST SNP-RISK ALLELE"].str.split('-').str[1],
        "cancer_type": cancer_type
    })
    
    # Sort by chromosome and position
    bed_df = bed_df.sort_values(["chrom", "start", "end"])
    bed_df = bed_df.drop_duplicates()
    # Write to a BED file
    bed_df.to_csv(output_file, sep='\t', header=False, index=False, mode='a')

# Input files and cancer types
input_files = [
    ("/Users/noah/Downloads/renal.gwas_catalog.01_17_25.filtered.annotated.tsv", "renal"),
    ("/Users/noah/Downloads/breast.gwas_catalog.01_17_25.filtered.annotated.tsv", "breast"),
    ("/Users/noah/Downloads/colorectal.gwas_catalog.01_17_25.filtered.annotated.tsv", "colorectal"),
    ("/Users/noah/Downloads/lung.gwas_catalog.01_17_25.filtered.annotated.tsv", "lung"),
    ("/Users/noah/Downloads/prostate.gwas_catalog.01_17_25.filtered.annotated.tsv", "prostate"),
]

# Process each file
output_file = "gwas_variants.bed"
for input_file, cancer_type in input_files:
    create_bed_file(input_file, cancer_type, output_file)
    #print(f"Processed {input_file} and saved as {output_file}")
