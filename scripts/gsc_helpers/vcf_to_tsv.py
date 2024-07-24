import pandas as pd
import numpy as np

# Function to parse VCF and create DataFrame
def vcf_to_tsv(vcf_file,risk_snp_dict):
    # Initialize lists and dictionaries
    data = {}
    sample_ids = []
    variants = []

    # Open VCF file and read line by line
    with open(vcf_file, 'r') as f:
        for line in f:
            line = line.strip()
            
            # Skip header lines starting with ##
            if line.startswith('##'):
                continue
            
            # Process header line starting with #CHROM
            elif line.startswith('#CHROM'):
                fields = line.split('\t')
                sample_ids = fields[9:]  # Extract sample IDs from header
                continue
            
            # Process variant data lines
            fields = line.split('\t')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alt = fields[4]
            allele = risk_snp_dict[f"{chrom}:{pos}"]
            key = f"{chrom}:{pos}-{allele}"
            variants.append(key)

            # Process genotype information for each sample
            for i, sample_id in enumerate(sample_ids):
                genotype_info = fields[i + 9]  # Start index for sample genotype info
                if './.' in genotype_info:
                    value = np.nan
                else:
                    alleles = genotype_info.split(':')[0].split('/')
                    if allele == ref:
                        value = 2 - (int(alleles[0]) + int(alleles[1]))
                    elif allele == alt:
                        value = int(alleles[0]) + int(alleles[1])
                    else:
                        value = 0
                    
                    if '.' in alleles:
                        value = np.nan
                    #else:
                         
                        #count_allele = alleles.count(str(allele_type))
                        #value = count_allele
                        #value = sum(int(allele) for allele in alleles)

                if sample_id not in data:
                    data[sample_id] = {}

                data[sample_id][key] = value

    # Create DataFrame
    df = pd.DataFrame.from_dict(data, orient='index')
    df.index.name = 'sample_id'
    df.columns = variants

    return df

def tsv_to_dict(tsv_file):
    data_dict = {}
    with open(tsv_file, 'r') as f:
        for line in f:
            if line.startswith('#'):  # Skip header lines if any
                continue
            fields = line.strip().split('\t')
            key = f"{fields[0]}:{fields[2]}"
            value = fields[3]
            data_dict[key] = value
    return data_dict

# Example usage
vcf_file = 'PROFILE.noncoding.vcf'
bed_file = 'noncoding_vars.hg19.bed'
output_tsv = 'output.tsv'

risk_snp_dict = tsv_to_dict(bed_file)
df = vcf_to_tsv(vcf_file,risk_snp_dict)
#risk_snp_dict = tsv_to_dict(bed_file)
df.to_csv(output_tsv, sep='\t', na_rep='NaN')

