# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This code will gather germline and somatic convergence.


version 1.0

import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/infer_sex_ancestry/wdl/Coding_Convergence.wdl" as GSC_CC

# The workflow gathers germline-somatic convergence for all combinations.
workflow Convergence {
  input{
    Array[File] germline_noncoding_vcfs
    Array[File] somatic_noncoding_vcfs
    Array[File] somatic_noncoding_vcf_idxs
    Array[File] somatic_tsv_arr
    Array[String] patient_ids
    File germline_merged_vep_vcf
    File germline_somatic_table
    File noncoding_germline_somatic_table
    File somatic_noncoding_hotspots_bed
    String cancer_type
  }
  
  scatter(patient_index in range(length(patient_ids))){
    call GSC_CC.convergence {
      input:
    	id = patient_ids[patient_index],
        cancer_type = cancer_type,
        somatic_tsv = somatic_tsv_arr[patient_index],
        germline_merged_vep_vcf = germline_merged_vep_vcf,
        germline_somatic_table = germline_somatic_table
    }
    
    call Extract_Noncoding_Germline_SNPs {
      input:
        patient_id = patient_ids[patient_index],
        cancer_type = cancer_type,
        germline_noncoding_vcf = germline_noncoding_vcfs[patient_index],
        noncoding_germline_somatic_table = noncoding_germline_somatic_table
    }
    
    call Extract_Somatic_Noncoding_Mutations {
      input:
        somatic_noncoding_vcf = somatic_noncoding_vcfs[patient_index],
        somatic_noncoding_vcf_idx = somatic_noncoding_vcf_idxs[patient_index],
        somatic_noncoding_hotspots_bed = somatic_noncoding_hotspots_bed,
        patient_id = patient_ids[patient_index],
        cancer_type = cancer_type
    }
    
    call Concatnate_Strings {
      input:
        coding_str = convergence.convergence,
        germline_noncoding_str = Extract_Noncoding_Germline_SNPs.out,
        somatic_noncoding_str = Extract_Somatic_Noncoding_Mutations.out
    }
  }
  
  call Write_tsv {
    input:
      cancer_type = cancer_type,
      patients_info = Concatnate_Strings.out,
      germline_somatic_table = germline_somatic_table,
      noncoding_germline_somatic_table = noncoding_germline_somatic_table,
      somatic_noncoding_hotspots_bed = somatic_noncoding_hotspots_bed
  }
  
  output{
    File out = Write_tsv.out
  }
}


# Task to Extract Noncoding Germline SNPs
task Extract_Noncoding_Germline_SNPs {
  input {
    File germline_noncoding_vcf  # Germline VCF file
    File noncoding_germline_somatic_table  # Table defining noncoding germline somatic regions
    String cancer_type  # Type of cancer being analyzed
    String patient_id  # Unique identifier for the sample
  }
  
  command <<<
    # Create an empty file to store extracted SNPs
    touch ~{patient_id}.vars.txt
    
    # Extract SNP IDs from the noncoding germline somatic table for the specified cancer type
    cut -f1,2 < ~{noncoding_germline_somatic_table} | tail -n +2 | grep -i '~{cancer_type}' | cut -f2 | sort | uniq > germline_noncoding.list
    
    # Initialize an empty string to store SNP information
    out=""
    
    # Extract CHROM, POS, REF, ALT, and GT fields from the VCF file
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' ~{germline_noncoding_vcf} | awk 'BEGIN {FS="\t"; OFS="\t"} {
    split($5, a, "/")
    if (a[1] == "1" && a[2] == "1") {
        print $1":"$2"-"$4, "HOM"
    } else if (a[1] == "0" && a[2] == "1") {
        print $1":"$2"-"$3, "HET"
        print $1":"$2"-"$4, "HET"
    } else if (a[1] == "1" && a[2] == "0") {
        print $1":"$2"-"$3, "HET"
        print $1":"$2"-"$4, "HET"
    } else if (a[1] == "0" && a[2] == "0") {
        print $1":"$2"-"$3, "HOM"
    }
    }' > query.txt
    
    # Filter the SNP table to get CHR, POS, RISK_ALLELE, then search through query.txt for it
    grep -Fwf germline_noncoding.list query.txt > ~{patient_id}.vars.txt
    
    # Loop through each SNP in the SNP list
    while IFS= read -r snp; do
      # Check if SNP is homozygous for risk allele
      if [ $(grep "$snp	HOM" ~{patient_id}.vars.txt | awk 'END {print NR}') -gt 0 ]; then
        out="${out}	2"  # Append '2' to the output string
      # Check if SNP is heterozygous
      elif [ $(grep "$snp	HET" ~{patient_id}.vars.txt | awk 'END {print NR}') -gt 0 ]; then
        out="${out}	1"  # Append '1' to the output string
      else
        out="${out}	0"  # Append '0' to the output string
      fi
    done < germline_noncoding.list

    # Write the output string to a file
    echo -e "$out" > ~{patient_id}.txt
  >>>
  
  output {
    File out = "~{patient_id}.txt"
  }
  runtime {
    docker: "vanallenlab/bcftools"
  }
}

# Task to Extract Somatic Mutations in Noncoding Regions
task Extract_Somatic_Noncoding_Mutations {
  input {
    File somatic_noncoding_vcf # Somatic VCF file
    File somatic_noncoding_vcf_idx # Index file for the somatic VCF
    File somatic_noncoding_hotspots_bed # BED file defining noncoding regions
    String patient_id  # Unique identifier for the patient (tumor(s) for each patient are previously put into one VCF)
    String cancer_type # Type of cancer the patient has
    
  }
  command <<<
    # Filter noncoding regions BED file to include only regions relevant to the specified cancer type
    grep -i '~{cancer_type}' ~{somatic_noncoding_hotspots_bed} | sort > somatic_noncoding_hotspots.bed
    
    # Initialize an empty string to store somatic mutation counts
    somatic_noncoding_str=""

    # Iterate through each noncoding regions
    while IFS=$'\t' read -r CHROM START END TISSUE GENE; do

      # Count somatic mutations within the region using bcftools
      somatic_mutation_count=$(bcftools view -H -r $CHROM:$START-$END ~{somatic_noncoding_vcf} | wc -l)

      # Append the mutation count to the somatic mutations string
      somatic_noncoding_str="$somatic_noncoding_str	$somatic_mutation_count"
    done < somatic_noncoding_hotspots.bed

    # Write somatic mutation counts to a text file
    echo "$somatic_noncoding_str" > noncoding_somatic_mutation_counts.txt
  >>>

  output {
    File out = "noncoding_somatic_mutation_counts.txt"
  }
  runtime {
    docker: "vanallenlab/bcftools"
  }
}

# Task to Extract Somatic Mutations in Noncoding Regions
task Write_tsv{
  input{
    Array[File] patients_info
    File germline_somatic_table
    File noncoding_germline_somatic_table
    File somatic_noncoding_hotspots_bed
    String cancer_type
  }
  
  command <<<
    #Only interested in coding variants
    cut -f1,2,3 < ~{germline_somatic_table} | tail -n +2 | grep -v 'noncoding' | grep -i '~{cancer_type}' | cut -f2 | sort | uniq > germline_coding.list
    cut -f1,7,8 < ~{germline_somatic_table} | tail -n +2 | grep -v 'noncoding' | grep -i '~{cancer_type}' | cut -f2 | sort | uniq > somatic_coding.list
    cut -f1,2 < ~{noncoding_germline_somatic_table} | tail -n +2 | grep -i '~{cancer_type}' | cut -f2 | sort | uniq > germline_noncoding.list
    grep -i '~{cancer_type}' ~{somatic_noncoding_hotspots_bed} | sort | cut -f5 > noncoding_somatic_regions.list
    
    header="patient_id"
    while IFS= read -r gene; do
      header="${header}	$gene-g"
    done < germline_coding.list
    while IFS= read -r gene; do
      header="${header}	$gene-s"
    done < somatic_coding.list
    while IFS= read -r snp; do
      header="${header}	$snp"
    done < germline_noncoding.list
    while IFS= read -r gene_region; do
      header="${header}	$gene_region-ncs"
    done < noncoding_somatic_regions.list
    
    echo "$header" > ~{cancer_type}.tsv
    cat ~{write_lines(patients_info)} >> patient_files.list
    while IFS= read -r file;do
      cat "$file" >> ~{cancer_type}.tsv
    done < patient_files.list
  >>>
  runtime {
    docker: "ubuntu:latest"
  }
  output {
  	File out = "~{cancer_type}.tsv"
  }
}

task Concatnate_Strings {
  input{
    File coding_str
    File germline_noncoding_str
    File somatic_noncoding_str
  }
  command <<<
    str1=$(cat ~{coding_str})
    str2=$(cat ~{germline_noncoding_str})
    str3=$(cat ~{somatic_noncoding_str})
    echo "$str1$str2$str3" > out.txt
  >>>
  runtime{
    docker: "ubuntu:latest"
  }
  output {
    File out = "out.txt"
  }
}

