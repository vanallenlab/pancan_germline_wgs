version 1.0

task Noncoding_germline_SNPs {
  input {
    File vcfFile
    String cancer_type
    String id
    File noncoding_germline_somatic_table
    
  }
  command <<<
    touch ~{id}.vars.txt
    cut -f1,2 < ~{noncoding_germline_somatic_table} | tail -n + 2 | grep -i '~{cancer_type}' | cut -f2 | sort | uniq > germline_noncoding.list
    
    out="~{id}"
    
    #Extract CHROM, POS, REF,ALT, and GT from vcf from
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT]\n' ~{vcfFile} | awk 'BEGIN {FS="\t"; OFS="\t"} {
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
    
    #Filter the bed file to get just CHR,POS,RISK_ALLELE; then search through query.txt for it.
    grep -Fwf germline_noncoding.list query.txt > ~{id}.vars.txt
    
    while IFS= read -r snp; do
      if [ $(grep "$snp	HOM" ~{id}.vars.txt) -gt 0 ];then
        out="{out}	2"
      elif [ $(grep "$snp	HET" ~{id}.vars.txt) -gt 0 ];then
        out="{out}	1"
      else
        out="{out}	0"
      fi
    done < germline_noncoding.list
    echo "$out" > ~{id}.snps
  >>>
  output {
    File out = "~{id}.snps"
  }
  runtime {
    docker: "vanallenlab/bcftools"
  }
}

workflow noncoding_portion {
  input {
    File germline_noncoding_vcfFile
    String id
    String cancer_type
    File noncoding_germline_somatic_table
  }
  
  call Noncoding_germline_SNPs {
    input:
      vcfFile = germline_noncoding_vcfFile,
      cancer_type = cancer_type,
      id = id,
      noncoding_germline_somatic_table = noncoding_germline_somatic_table
  }
  output {
    File noncoding_patient_info = Noncoding_germline_SNPs.out
  }
}
