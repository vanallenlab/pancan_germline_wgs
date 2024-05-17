#!/bin/bash


#Check if a vcf file is sorted

# Define function to check if variants are properly sorted
check_sorted() {
    local vcf_file="$1"

    # Extract chromosome and position from VCF file
    bcftools query -f '%CHROM\t%POS\n' "${vcf_file}" > chrom_pos.txt

    # Check if variants are properly sorted
    if sort -c chrom_pos.txt; then
        echo "VCF file is properly sorted."
    else
        echo "VCF file is not properly sorted."
    fi

    # Clean up temporary file
    #rm chrom_pos.txt
}

# Main script
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <vcf_file>"
    exit 1
fi

vcf_file="$1"
check_sorted "${vcf_file}"