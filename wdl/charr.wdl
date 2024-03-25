# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fields@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# The purpose of this code is to use CHARR, a 
# contamination estimator which leverages the 
# infiltration of reference reads within homozygous 
# alternate variant calls.
# https://github.com/atgu/CHARR
# We also use echtvar for annotation of gVCF
# https://github.com/brentp/echtvar


version 1.0

task decompose {
  input {
  	File reblocked_gvcf #g.vcf.gz file
    Int? disk_gb
    Float memory_var = 3.5
  }
  Int default_disk_size = ceil((2 * size(reblocked_gvcf, "GB")) + 10)
  
  command <<<
    bcftools norm -m - ~{reblocked_gvcf} -o clean_bcf.bcf
  >>>
  
  runtime {
    docker: "vanallenlab/charr"
    memory: "~{memory_var} GB"  # Adjust the amount based on your requirements
    disks: "local-disk " + select_first([disk_gb, default_disk_size]) + " HDD"
  }
  output {
  	File clean_bcf = "clean_bcf.bcf"
  }
}

#We use echtvar to annotate
task annotate {
  input {
  	File clean_bcf
    File gnomad_ref
    Int? disk_gb
    Float memory_var = 3.5
  }
  Int default_disk_size = ceil((2 * size(clean_bcf, "GB")) + 10)
  
  command <<<
    echtvar anno -e ~{gnomad_ref} ~{clean_bcf} output.annotated.bcf
  >>>  
  
  runtime {
    docker: "vanallenlab/charr"
    memory: "~{memory_var} GB"  # Adjust the amount based on your requirements
    disks: "local-disk " + select_first([disk_gb, default_disk_size]) + " HDD"
  }
  output {
  	File annotated_bcf = "output.annotated.bcf"
  }
}
task step1 {
  input {
  	File gvcf #g.vcf.gz file
    Int? disk_gb
    Float memory_var = 3.5
  }
  Int default_disk_size = ceil((2 * size(gvcf, "GB")) + 10)
  
  command <<<
    bcftools view -i 'ALT!="<NON_REF>" & SUM(FORMAT/AD) > 0' ~{gvcf} > step1.annotated.bcf
  >>>
  
  runtime {
    docker: "vanallenlab/charr"
    memory: "~{memory_var} GB"  # Adjust the amount based on your requirements
    disks: "local-disk " + select_first([disk_gb, default_disk_size]) + " HDD"
  }
  output {
  	File output_bcf = "step1.annotated.bcf"
  }
}

task step2 {
  input {
  	File gvcf #g.vcf.gz file
    Int? disk_gb
    Float memory_var = 3.5
  }
  Int default_disk_size = ceil((2 * size(gvcf, "GB")) + 10)
  
  
  command <<<
    bcftools annotate -c INFO/AF:=INFO/gnomad_popmax_af ~{gvcf} > step2.annotated.bcf
  >>>
  
  runtime {
    docker: "vanallenlab/charr"
    memory: "~{memory_var} GB"  # Adjust the amount based on your requirements
    disks: "local-disk " + select_first([disk_gb, default_disk_size]) + " HDD"
  }
  output {
  	File output_bcf = "step2.annotated.bcf"
  }
}

task contamination {
  input {
  	File gvcf #g.vcf.gz file
    Int? disk_gb
    Float memory_var = 3.5
  }
  
  Int default_disk_size = ceil((2 * size(gvcf, "GB")) + 10)
  
  command <<<
    sceVCF -o contamination.txt ~{gvcf}
  >>>
  
  runtime {
    docker: "vanallenlab/charr"
    memory: "~{memory_var} GB"  # Adjust the amount based on your requirements
    disks: "local-disk " + select_first([disk_gb, default_disk_size]) + " HDD"
  }
  output {
  	File charr_output = "contamination.txt"
  }
}

# Define workflow
workflow charr {
  input {
  	 File gnomad_ref = "gs://dfci-g2c-refs/echtvar/gnomad.v3.1.2.echtvar.popmax.v2.zip"
     File reblocked_gvcf
  }
 
  
  call decompose{
  	input:
      reblocked_gvcf = reblocked_gvcf
  }
  call annotate{
    input: 
      clean_bcf = decompose.clean_bcf,
      gnomad_ref = gnomad_ref
  }
  call step1{
  	input:
    	gvcf = annotate.annotated_bcf
  }
  call step2{
  	input:
    	gvcf = step1.output_bcf
  }
  call contamination{
    input:
      gvcf = step2.output_bcf
  }
  
  output{
    File charr_output = contamination.charr_output
  }
}