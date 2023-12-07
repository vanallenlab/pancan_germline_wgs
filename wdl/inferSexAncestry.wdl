# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

#The purpose of this code is to infer sex and ancestry
#from variant call data. 


version 1.0

task get_chr_count {
  input {
  	File counts_coverage_file
  }

  command <<<
    Rscript /opt/pancan_germline_wgs/scripts/ufc_helpers/sex_infer.R ~{counts_coverage_file} > log.txt
    
    head -n1 log.txt | cut -f1 -d'-' > chrX.txt
    head -n1 log.txt | cut -f2 -d'-' > chrY.txt

    #echo "Result: ${result}" > log.txt
    
    # Extract numbers from the result
    #echo 4 > chrX.txt
    #echo 5 > chrY.txt
    #echo "${result}" | cut -d'-' -f1 > chrX.txt
    #echo "${result}" | cut -d'-' -f2 > chrY.txt

    # Write the numbers to respective files
    #echo "$chrX" > chrX.txt
    #echo "$chrY" > chrY.txt
  >>>
  
  runtime {
    docker: "vanallenlab/g2c_pipeline:nf_sex_ancestry_alpha"
    memory: "8 GB"  # Adjust the amount based on your requirements
  }
  output {
    Int chrX = read_int("chrX.txt")
    Int chrY = read_int("chrY.txt")
    File output_log = "log.txt"
  }
}

task infer_ancestry{
  input {
    File gvcf_file
  }
  command <<<
    cp /AncInferSNPs.txt .
    grafpop ~{gvcf_file} output.txt
    /opt/infer_ancestry.py output.txt > tmp.txt #outputs to ancestryID.txt and ancestry.txt
    
    # Extract the first line and save it to 'ancestry.txt'
    head -n 1 tmp.txt > ancestry.txt

    # Extract the second line and save it to 'ancestryID.txt'
    tail -n 1 tmp.txt > ancestryID.txt

    
    
  >>>
  runtime{
    docker: "vanallenlab/g2c_pipeline:nf_graf_plink_alpha"
    memory: "60 GB"  # Adjust the amount based on your requirements
    bootDiskSizeGb: 10 # Adjust the amount based on your requirements
  }
  output{
    Int Ancestry_ID = read_int("ancestryID.txt")
    String Ancestry = read_string("ancestry.txt")

    
  }
}

# Define workflow
workflow infer_sex {

  input {
    File counts_coverage_file
    File gvcf_file
    
  }
  
  call get_chr_count{
    input:
      counts_coverage_file = counts_coverage_file,
  }
  call infer_ancestry{
    input:
      gvcf_file = gvcf_file
  }
  output{
  	Int chrX_count = get_chr_count.chrX
    Int chrY_count = get_chr_count.chrY
    Int Ancestry_ID = infer_ancestry.Ancestry_ID
    String Ancestry = infer_ancestry.Ancestry
    #File ls2_log = infer_ancestry.ls2_log
    #File ls_log = infer_ancestry.ls_log
    
  }
}

