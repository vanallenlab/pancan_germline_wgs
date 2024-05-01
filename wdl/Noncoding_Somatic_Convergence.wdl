version 1.0

task Noncoding_Somatic_Regions {
  input {
    File somatic_noncoding_vcfFile
    File noncoding_regions_bed
    String id
    String cancer_type
    
  }
  command <<<
    grep -i '~{cancer_type}' ~{noncoding_regions_bed} | sort > noncoding_somatic_regions.bed
    
    somatic_nc_str="~{id}"
    while IFS= read -r CHROM START END TISSUE GENE; do
      som_mut_count=$(bcftools view -H -r "$CHROM:$START-$END" ~{somatic_noncoding_vcfFile} | wc -l)
      somatic_nc_str="$somatic_nc_str\t$som_mut_count"
    done < noncoding_somatic_regions.bed
    echo "$somatic_nc_str" > somatic_mutations.txt
  >>>

  output {
    File out = "somatic_mutations.txt"
  }
  runtime {
    docker: "vanallenlab/bcftools"
  }
}

workflow gather_nc_somatic_mutations {
  input {
    File somatic_noncoding_vcfFile
    String id
    String cancer_type
    File noncoding_regions_bed
  }
  
  call Noncoding_Somatic_Regions {
    input:
      somatic_noncoding_vcfFile = somatic_noncoding_vcfFile,
      cancer_type = cancer_type,
      id = id,
      noncoding_regions_bed = noncoding_regions_bed
  }
  output {
    File somatic_noncoding_info = Noncoding_Somatic_Regions.out
  }
}
