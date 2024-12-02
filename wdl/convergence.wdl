# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

version 1.0

# Analysis of germline variants and somatic drivers in coding regions
# This workflow identifies pathogenic and high-impact variants in germline data, as well as somatic driver mutations in coding regions.

# Task to extract coding regions from a somatic table
task Extract_Coding_Regions {
  input {
    File germline_somatic_table
    String cancer_type
  }
  command <<<
    #Only interested in coding variants
    cut -f1,2,3,4 < ~{germline_somatic_table} | tail -n +2 | grep -v 'noncoding' | grep -v 'UNK' | grep -i '~{cancer_type}' | cut -f2 | sort | uniq > germline_coding.list
    cut -f1,7,8 < ~{germline_somatic_table} | tail -n +2 | grep -v 'noncoding' | grep -i '~{cancer_type}' | cut -f2 | sort | uniq > somatic_coding.list
  >>>
  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }
  output {
  	File germline_coding_list = "germline_coding.list"
    File somatic_coding_list = "somatic_coding.list"
  }
}

# Task to extract pathogenic and high-impact germline variants
task Extract_Germline_Variants {
  input {
    File germline_merged_vep_vcf
    File germline_somatic_table
    File germline_genes
    String patient_id
  }
  
  command <<<
  set -eu -o pipefail

  # Extract variants present only in the germline data
  bcftools view -s ~{patient_id} ~{germline_merged_vep_vcf} -o sample.vcf
  bcftools view -i 'AC>0 & GT="alt"' sample.vcf -o germline_only.vcf
  echo "Checkpoint 1"

  # Split VEP annotations and filter for relevant information
	bcftools +split-vep -f '%CHROM|%POS|%CSQ/SYMBOL|%CSQ/IMPACT|%CSQ/gnomAD_AF_nfe|%CSQ/gnomAD_AF_popmax|%CSQ/gnomAD_controls_AF_popmax|%CSQ/ClinVar_external_CLNSIG|%CSQ/Consequence' germline_only.vcf | cut -d'|' -f1,2,5,6,39-44 | awk '{gsub(/\|/, "\t");print}' | cut -d',' -f1 | awk -F'\t' '($5 == "Pathogenic" || $5 == "Likely_pathogenic" || $9 < 0.02 || $9 == "." || $9 == "") &&($5 == "Pathogenic" || $5 == "Likely_pathogenic" || $7 < 0.02 || $7 == "." || $7 == "") && ($5 == "Pathogenic" || $5 == "Likely_pathogenic" || $8 < 0.02 || $8 == "." || $8 == "")' | grep -Ev 'Benign|Likely_benign|LOW' | grep -Fwf ~{germline_genes} > query.tsv
  
  echo "Checkpoint 2"
  # Extract potential deleterious germline genes
  cut -f4 query.tsv | sort | uniq > potential_deleterious_germline_genes.list
  touch c1.list
  touch c2.list
  touch c3.list
  touch pathogenic_germline_genes.list

  python3 <<CODE
  import pandas as pd

  # Load data
  vep_df = pd.read_csv("query.tsv", sep='|', names=["CHROM", "POS", "Gene", "IMPACT", "gnomad_AF_nfe", "gnomad_AF_popmax", "gnomAD_controls_AF_popmax", "ClinVar_CLNSIG", "Consequence"])
  convergence_table = pd.read_csv("~{germline_somatic_table}", sep='\t')

  # Create arrays for gene categories
  tsg_arr = list(set(convergence_table[convergence_table['germline_onc_status'] == "TSG"]['germline_gene'].tolist()))
  onc_tsg_arr = list(set(convergence_table[convergence_table['germline_onc_status'] == "ONC-TSG"]['germline_gene'].tolist()))

  # Create arrays for mutation types
  deletion_arr = list(set(convergence_table[convergence_table['germline_mutations'].str.contains('D', na=False)]['germline_gene'].tolist()))
  nonsense_arr = list(set(convergence_table[convergence_table['germline_mutations'].str.contains('N', na=False)]['germline_gene'].tolist()))
  missense_arr = list(set(convergence_table[convergence_table['germline_mutations'].str.contains('Mis', na=False)]['germline_gene'].tolist()))
  splice_site_arr = list(set(convergence_table[convergence_table['germline_mutations'].str.contains('S', na=False)]['germline_gene'].tolist()))
  translocation_arr = list(set(convergence_table[convergence_table['germline_mutations'].str.contains('T', na=False)]['germline_gene'].tolist()))
  frameshift_arr = list(set(convergence_table[convergence_table['germline_mutations'].str.contains('F', na=False)]['germline_gene'].tolist()))

  # Initialize an empty set for pathogenic genes
  pathogenic_gene_list = set()
  print("Checkpoint 3")
  
  # Loop through each row in the DataFrame
  for _, row in vep_df.iterrows():
      # If there is a P/LP variant, add the gene
      if row['ClinVar_CLNSIG'] in ['Pathogenic', "Likely_pathogenic", "Pathogenic/Likely_pathogenic"]:
          pathogenic_gene_list.add(row['Gene'])
          continue

      # Check for HIGH-impact mutations in TSG or ONC-TSG genes
      if (row['Gene'] in tsg_arr or row['Gene'] in onc_tsg_arr) and \
         (row['Consequence'].lower() in ["deletion", "frameshift", "nonsense"] and row['IMPACT'] == "HIGH"):
          pathogenic_gene_list.add(row['Gene'])
          continue

  # Write the list to a file
  with open("~{patient_id}.txt", "w") as file:
      for gene in sorted(pathogenic_gene_list):  # Sorting is optional but improves readability
          file.write(f"{gene}\n")

  print(f"Pathogenic genes written to ~{patient_id}.txt")
  CODE

  >>>
  
  runtime {
    docker: "vanallenlab/g2c_pipeline"
    preemptible: 3
  }
  output{
  	File out = "~{patient_id}.txt"
  }
}

# Extract somatic driver mutations based on somatic genes
task Extract_Somatic_Drivers {
  input {
	File somatic_tsv
	File somatic_genes
  String patient_id
  }
  command <<<
	# Extract somatic driver mutations based on somatic genes
	cat ~{somatic_tsv} | cut -f3 | grep -Fwf ~{somatic_genes} | sort | uniq >> ~{patient_id}.somatic_drivers.txt
  >>>
  
  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }
  output{
  	File out = "~{patient_id}.somatic_drivers.txt"
  }

}

task Output_Mutations {
  input {
    File germline_variants
    File somatic_drivers
    String patient_id
    File germline_genes
    File somatic_genes
  }
  
  command <<<
	# Initialize output string with patient ID
	output="~{patient_id}"

	# Check each germline gene for pathogenic variants
	while IFS=$'\t' read -r germline_gene; do
		if [ $(grep -F "$germline_gene" ~{germline_variants} | awk 'END {print NR}' ) -gt 0 ]; then
			output="${output}	1"
		else
			output="${output}	0"
		fi
	done < ~{germline_genes}

	# Check each somatic gene for driver mutations
	while IFS=$'\t' read -r somatic_gene; do
		if [ $(grep -F "$somatic_gene" ~{somatic_drivers} | awk 'END {print NR}' ) -gt 0 ]; then
			output="${output}	1"
		else
			output="${output}	0"
		fi
	done < ~{somatic_genes}

	# Write output string to file
	echo "$output" > ~{patient_id}_info.txt
  >>>
  
  runtime {
    docker: "ubuntu:latest"
    preemptible: 3
  }
  output{
  	File out = "~{patient_id}_info.txt"
  }
}


# Define workflow
workflow convergence {
  input{
    File germline_merged_vep_vcf # Merged VCF file annotated with VEP
    File germline_somatic_table # Germline-Somatic table file containing convergence pairs
    File somatic_tsv # Somatic TSV file containing mutation information
    String patient_id # Unique identifier for the patient
    String cancer_type # Type of cancer being analyzed
  }
  call Extract_Coding_Regions {
    input:
      germline_somatic_table = germline_somatic_table,
      cancer_type = cancer_type
  }
  call Extract_Germline_Variants {
    input:
      germline_merged_vep_vcf = germline_merged_vep_vcf,
      germline_genes = Extract_Coding_Regions.germline_coding_list,
      germline_somatic_table = germline_somatic_table,
      patient_id = patient_id
  }
  call Extract_Somatic_Drivers {
    input:
      somatic_genes = Extract_Coding_Regions.somatic_coding_list,
      somatic_tsv = somatic_tsv,
      patient_id = patient_id
  }
  call Output_Mutations {
	input:
      patient_id = patient_id,
      somatic_drivers = Extract_Somatic_Drivers.out,
      germline_variants = Extract_Germline_Variants.out,
      germline_genes = Extract_Coding_Regions.germline_coding_list,
	    somatic_genes = Extract_Coding_Regions.somatic_coding_list
  }
  output{
  	File convergence = Output_Mutations.out
  }
}