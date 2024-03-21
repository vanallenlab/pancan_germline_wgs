# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# This code will standardize a HMF vcf.


version 1.0

task standardize{
  input {
	File unstandardized_vcf
	String id
  }
  Int default_disk_size = ceil((2 * size(unstandardized_vcf, "GB")) + 10)
  
  command <<<
	# Step 1: Remove annotations (INFO fields)
	bcftools view -h ~{unstandardized_vcf} | tail -n 1 | cut -f10 > sample_name.txt
	bcftools annotate -x FORMAT/AD,INFO ~{unstandardized_vcf} -o temp.vcf.gz

	# Step 2: Keep only samples listed in sample_name.txt
	bcftools view -S sample_name.txt temp.vcf.gz -o filtered_samples.vcf.gz

	# Step 3: Normalize the VCF file
	bcftools norm -m - filtered_samples.vcf.gz -o normalized.vcf.gz

	# Step 4: Filter variants where 'AC > 0' and the filter passes QC
	bcftools filter -i 'AC > 0 & FILTER == "PASS"' normalized.vcf.gz -o allele_count.vcf.gz

	# Step 5:
	bcftools +fill-tags allele_count.vcf.gz -o ~{id}.vcf.gz -- -t FORMAT/VAF

	# Step 6: Index the vcf
	bcftools index -t -o ~{id}.vcf.gz.tbi ~{id}.vcf.gz
  >>>
  
  runtime {
    docker: "vanallenlab/bcftools"
    disks: "local-disk " + default_disk_size + " HDD"
  }
  output{
    File vcf = "~{id}.vcf.gz"
    File vcf_idx  = "~{id}.vcf.gz.tbi"
  }
}




# Define workflow
workflow standardize_hmf {
  input{
	File unstandardized_vcf
    String id
  }
  call standardize {
	input:
		id = id,
		unstandardized_vcf = unstandardized_vcf
  }
  output{
	File standardized_vcf = standardize.vcf
    File standardized_vcf_idx = standardize.vcf_idx
  }
}