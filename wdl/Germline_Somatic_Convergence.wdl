# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This code will gather germline and somatic convergence.


version 1.0

task get_coding_regions{
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

task get_germline_coding_variants{
  input {
    File germline_merged_vep_vcf
    File germline_somatic_table
    File germline_genes
    String id
  }
  
  command <<<
	bcftools view -s ~{id} ~{germline_merged_vep_vcf} -o sample.vcf
	bcftools view -i 'AC>0 & GT="alt"' sample.vcf -o germline_only.vcf

	bcftools +split-vep -f '%CHROM|%POS|%CSQ/SYMBOL|%CSQ/IMPACT|%CSQ/gnomAD_AF_nfe|%CSQ/gnomAD_AF_popmax|%CSQ/gnomAD_controls_AF_popmax|%CSQ/ClinVar_external_CLNSIG' germline_only.vcf | cut -d'|' -f1,2,5,6,39-43 | awk '{gsub(/\|/, "\t"); print}' | cut -d',' -f1 | awk -F'\t' '($9 < 0.02 || $9 == "." || $9 == "") && ($7 < 0.02 || $7 == "." || $7 == "") && ($8 < 0.02 || $8 == "." || $8 == "")' | grep -Fwf ~{germline_genes} > query.tsv
	cut -f4 query.tsv | sort | uniq > potential_deleterious_germline_genes.list
	touch c1.list
	touch c2.list
	touch c3.list
	touch pathogenic_germline_genes.list

	# Condition 1: Denoted P/LP by ClinVar
	grep -E 'Pathogenic|Likely_pathogenic' query.tsv | cut -f4 | sort | uniq > c1.list

	# Condition 2: If includes TSG or Mutation-Type includes D,N,F
	while IFS= read -r gene; do
		TSG=$(cut -f2-5 ~{germline_somatic_table} | grep "$gene" | grep 'TSG' | awk 'END {print NR}')
		MT=$(cut -f2-5 ~{germline_somatic_table} | grep "$gene" | cut -f5 | grep -E 'D|N|F' | awk 'END {print NR}')
		if ((TSG + MT > 0)); then
			echo "$gene" >> c2.list
		fi
	done < potential_deleterious_germline_genes.list

	# Condition 3: If includes ONC or Mutation-Type includes Mis
	while IFS= read -r gene; do
		ONC=$(cut -f2-5 ~{germline_somatic_table} | grep "$gene" | grep 'ONC' | awk 'END {print NR}')
		MIS=$(cut -f2-5 ~{germline_somatic_table} | grep "$gene" | cut -f5 | grep 'Mis' | awk 'END {print NR}')
		if ((ONC + MIS > 0)); then
			echo "$gene" >> c3.list
		fi
	done < potential_deleterious_germline_genes.list

	# Deciding which genes are pathogenic
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

	sort pathogenic_germline_genes.list | uniq > ~{id}.txt


  >>>
  
  runtime {
    docker: "vanallenlab/bcftools"
  }
  output{
  	File variants = "~{id}.txt"
  }
}

task get_somatic_coding_drivers{
  input {
	File somatic_tsv
	File somatic_genes
    String id
  }
  command <<<
	cat ~{somatic_tsv} | cut -f3 | grep -Fwf ~{somatic_genes} | sort | uniq >> ~{id}.som_driving.txt
  >>>
  
  runtime {
    docker: "ubuntu:latest"
  }
  output{
  	File drivers = "~{id}.som_driving.txt"
  }

}

task converge{
  input {
    File germline_variants
    File somatic_drivers
    String id
    File germline_genes
    File somatic_genes
  }
  
  command <<<
	output="~{id}"
	while IFS=$'\t' read -r germline_gene; do
		if [ $(grep -F "$germline_gene" ~{germline_variants} | awk 'END {print NR}' ) -gt 0 ]; then
			output="${output}	1"
		else
			output="${output}	0"
		fi
	done < ~{germline_genes}

	while IFS=$'\t' read -r somatic_gene; do
		if [ $(grep -F "$somatic_gene" ~{somatic_drivers} | awk 'END {print NR}' ) -gt 0 ]; then
			output="${output}	1"
		else
			output="${output}	0"
		fi
	done < ~{somatic_genes}
	echo "$output" > ~{id}_info.txt
  >>>
  
  runtime {
    docker: "ubuntu:latest"
  }
  output{
  	File patient_info = "~{id}_info.txt"
  }
}


# Define workflow
workflow convergence {
  input{
    File germline_merged_vep_vcf
    File germline_somatic_table
    String id
    String cancer_type
    File somatic_tsv
  }
  call get_coding_regions{
    input:
      germline_somatic_table = germline_somatic_table,
      cancer_type = cancer_type
  }
  call get_germline_coding_variants{
    input:
      germline_merged_vep_vcf = germline_merged_vep_vcf,
      germline_genes = get_coding_regions.germline_coding_list,
      germline_somatic_table = germline_somatic_table,
      id = id
  }
  call get_somatic_coding_drivers{
    input:
      somatic_genes = get_coding_regions.somatic_coding_list,
      somatic_tsv = somatic_tsv,
      id = id
  }
  call converge{
	input:
      id = id,
      somatic_drivers = get_somatic_coding_drivers.drivers,
      germline_variants = get_germline_coding_variants.variants,
      germline_genes = get_coding_regions.germline_coding_list,
	  somatic_genes = get_coding_regions.somatic_coding_list
  }
  output{
  	File convergence = converge.patient_info
  }
}