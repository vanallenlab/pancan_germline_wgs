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
    cut -f1,2,3 < ~{germline_somatic_table} | tail -n +2 | grep -v 'noncoding' | grep -i '~{cancer_type}' | cut -f2 | sort | uniq > germline_coding.list
    cut -f1,7,8 < ~{germline_somatic_table} | tail -n +2 | grep -v 'noncoding' | grep -i '~{cancer_type}' | cut -f2 | sort | uniq > somatic_coding.list
  >>>
  runtime {
    docker: "ubuntu:latest"
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
	# Extract variants present only in the germline data
	bcftools view -s ~{patient_id} ~{germline_merged_vep_vcf} -o sample.vcf
	bcftools view -i 'AC>0 & GT="alt"' sample.vcf -o germline_only.vcf

	# Split VEP annotations and filter for relevant information
	bcftools +split-vep germline_only.vcf -f '%CHROM\t%POS\t%SYMBOL \
	\t%IMPACT\t%gnomAD_AF_nfe\t%gnomAD_AF_popmax\t%gnomAD_controls_AF_popmax\t%ClinVar_external_CLNSIG' -d \
	| awk -F'\t' '($8 == "Pathogenic" || $8 == "Likely_pathogenic" || $5 < 0.02 || $5 == "." || $5 == "") && ($8 == "Pathogenic" || $8 == "Likely_pathogenic" || $6 < 0.02 || $6 == "." || $6 == "") && ($8 == "Pathogenic" || $8 == "Likely_pathogenic" || $7 < 0.02 || $7 == "." || $7 == "")' | \ 
	grep -Ev 'Benign|Likely_benign|MODIFIER|LOW' | grep -Fwf ~{germline_genes} > query.tsv

	# Extract potential deleterious germline genes
	cut -f4 query.tsv | sort | uniq > potential_deleterious_germline_genes.list
	touch c1.list
	touch c2.list
	touch c3.list
	touch pathogenic_germline_genes.list

	# Condition 1: Pathogenic or Likely pathogenic by ClinVar
	grep -E 'Pathogenic|Likely_pathogenic' query.tsv | cut -f4 | sort | uniq > c1.list

	# Condition 2: Tumor suppressor gene (TSG) or Mutation-Type includes D,N,F
	while IFS= read -r gene; do
		TSG=$(cut -f2-5 ~{germline_somatic_table} | grep "$gene" | grep 'TSG' | awk 'END {print NR}')
		MT=$(cut -f2-5 ~{germline_somatic_table} | grep "$gene" | cut -f5 | grep -E 'D|N|F' | awk 'END {print NR}')
		if ((TSG + MT > 0)); then
			echo "$gene" >> c2.list
		fi
	done < potential_deleterious_germline_genes.list

	# Condition 3: Oncogene (ONC) or Mutation-Type includes Mis
	while IFS= read -r gene; do
		ONC=$(cut -f2-5 ~{germline_somatic_table} | grep "$gene" | grep 'ONC' | awk 'END {print NR}')
		MIS=$(cut -f2-5 ~{germline_somatic_table} | grep "$gene" | cut -f5 | grep 'Mis' | awk 'END {print NR}')
		if ((ONC + MIS > 0)); then
			echo "$gene" >> c3.list
		fi
	done < potential_deleterious_germline_genes.list

	# Determine pathogenic germline genes based on conditions
	while IFS= read -r gene; do
		# If condition 1 is met
		if [ $(grep "$gene" c1.list | awk 'END {print NR}') -gt 0 ];then
			echo "$gene" >> pathogenic_germline_genes.list
		# If condition 2 but not 3 is met
		elif [ $(grep "$gene" c2.list | awk 'END {print NR}') -gt 0 ] && [ $(grep "$gene" c3.list | awk 'END {print NR}') -eq 0 ]; then
			if [ $(grep "$gene" query.tsv | grep 'HIGH' | awk 'END {print NR}') -gt 0 ]; then
				echo "$gene" >> pathogenic_germline_genes.list
			fi
		# If condition 3 but not 2 is met
		elif [ $(grep "$gene" c2.list | awk 'END {print NR}') -eq 0 ] && [ $(grep "$gene" c3.list | awk 'END {print NR}') -gt 0 ]; then
			if [ $(grep "$gene" query.tsv | grep 'MODERATE' | awk 'END {print NR}') -gt 0 ]; then
				echo "$gene" >> pathogenic_germline_genes.list
			fi
		# If condition 2 and 3 is met or none of the conditions are met
		else
			if [ $(grep "$gene" query.tsv | grep -E 'HIGH|MODERATE' | awk 'END {print NR}') -gt 0 ]; then
				echo "$gene" >> pathogenic_germline_genes.list
			fi
		fi   	
	done < potential_deleterious_germline_genes.list

	sort pathogenic_germline_genes.list | uniq > ~{patient_id}.txt


  >>>
  
  runtime {
    docker: "vanallenlab/bcftools"
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