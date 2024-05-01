version 1.0

task Noncoding_Somatic_Regions {
  input {
    File somatic_noncoding_vcfFile
    File noncoding_regions_bed
    String id
    String cancer_type
    
  }
  command <<<
    grep -i '~{cancer_type}' ~{noncoding_regions_bed} > noncoding_somatic_regions.bed
    
    while IFS= read -r $CHROM $START $END $TISSUE $GENE; then
      som_mut=$(bcftools view -H -r "$CHROM:$START-$END" ~{somatic_noncoding_vcfFile} | wc -l)
      echo "$GENE\t$som_mut" >> somatic_mutations.txt
    done < noncoding_somatic_regions.bed

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
