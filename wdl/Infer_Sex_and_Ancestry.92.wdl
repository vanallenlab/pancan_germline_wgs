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

task organize {
  input {
  	File raw_ancestry_file
    String sample_ID
    Int ancestryID
    String ancestry
    Int chrX_count
    Int chrY_count
  }

  command <<<
    /opt/pancan_germline_wgs/scripts/ufc_helpers/organize_results.py --file ~{raw_ancestry_file} --ancestry ~{ancestry} --ancestryID ~{ancestryID} --chrX_count ~{chrX_count} --chrY_count ~{chrY_count} > ~{sample_ID}_demo.txt
  >>>
  
  runtime {
    docker: "vanallenlab/g2c_pipeline:nf_sex_ancestry_alpha"
    memory: "8 GB"  # Adjust the amount based on your requirements
  }
  output {
    File output_file = "~{sample_ID}_demo.txt"
  }
}

task infer_ancestry{
  input {
    File gvcf_file
    #Int memory_var = 60
    #Int? bootDiskSizeGb_var
    String sample_ID
  }
  Int memory_var = ceil(25 * size(gvcf_file, "GB"))
  Int bootDiskSizeGb_var = ceil(6 * size(gvcf_file, "GB")+10)
  
  command <<<
    cp /AncInferSNPs.txt .
    grafpop ~{gvcf_file} ~{sample_ID}_ancestry_raw.txt
    /opt/infer_ancestry.py ~{sample_ID}_ancestry_raw.txt > tmp.txt #outputs to ancestryID.txt and ancestry.txt
    
    # Extract the first line and save it to 'ancestry.txt'
    head -n 1 tmp.txt > ancestry.txt

    # Extract the second line and save it to 'ancestryID.txt'
    tail -n 1 tmp.txt > ancestryID.txt

  >>>
  runtime{
    docker: "vanallenlab/g2c_pipeline:nf_graf_plink_alpha"
    memory: "~{memory_var} GB"  # Adjust the amount based on your requirements
    bootDiskSizeGb: bootDiskSizeGb_var # Adjust the amount based on your requirements
  }
  output{
    Int Ancestry_ID = read_int("ancestryID.txt")
    String Ancestry = read_string("ancestry.txt")
    File grafpop_output_file = "~{sample_ID}_ancestry_raw.txt"

    
  }
}

# Define workflow
workflow infer_sex_ancestry {

  input {
    File counts_coverage_file
    File gvcf_file
    String sample_ID
    
  }
  
  call get_chr_count{
    input:
      counts_coverage_file = counts_coverage_file,
  }
  call infer_ancestry{
    input:
      gvcf_file = gvcf_file,
      sample_ID = sample_ID
  }
  call organize{
  	input:
      raw_ancestry_file = infer_ancestry.grafpop_output_file,
      sample_ID = sample_ID,
      ancestryID = infer_ancestry.Ancestry_ID,
      ancestry = infer_ancestry.Ancestry,
      chrX_count = get_chr_count.chrX,
      chrY_count = get_chr_count.chrY,
  }
  
  output{
    File demographics_file = organize.output_file
  }
}

