ref_table="../data/VALab_germline_somatic.tsv"


#Breast Section
grep 'breast' "$ref_table" | cut -f2,3 | grep 'noncoding' | cut -f1  > breast_noncoding_germline.tmp.list
cut -f3,10,15,18 breast_gwas_catalog_sig_filtered.annotated.tsv | tail -n +2 | grep -v 'coding_exon' | grep -Fwf breast_noncoding_germline.tmp.list | awk '{print $0 "\tbreast"}' > breast_noncoding_snps.tmp.tsv

#Colorectal Section
grep 'colorectal' "$ref_table" | cut -f2,3 | grep 'noncoding' | cut -f1 > colorectal_noncoding_germline.tmp.list
cut -f2,9,14,17 colorectal_gwas_catalog_sig_filtered.annotated.tsv | tail -n +2 | grep -v 'coding_exon' | grep -Fwf colorectal_noncoding_germline.tmp.list | awk '{print $0 "\tcolorectal"}' > colorectal_noncoding_snps.tmp.tsv

#Lung Section
grep 'lung' "$ref_table" | cut -f2,3 | grep 'noncoding' | cut -f1 > lung_noncoding_germline.tmp.list
cut -f2,9,14,17 lung_gwas_catalog_sig_filtered.annotated.tsv | tail -n +2 | grep -v 'coding_exon' | grep -Fwf lung_noncoding_germline.tmp.list | awk '{print $0 "\tlung"}' > lung_noncoding_snps.tmp.tsv

#Kidney Section
grep 'renal' "$ref_table" | cut -f2,3 | grep 'noncoding' | cut -f1 > kidney_noncoding_germline.tmp.list
cut -f2,9,14,17 renal_gwas_catalog_sig_filtered.annotated.tsv | tail -n +2 | grep -v 'coding_exon' | grep -Fwf kidney_noncoding_germline.tmp.list | awk '{print $0 "\trenal"}' > kidney_noncoding_snps.tmp.tsv

#Prostate Section
grep 'prostate' "$ref_table" | cut -f2,3 | grep 'noncoding' | cut -f1 > prostate_noncoding_germline.tmp.list
cut -f2,9,14,17 prostate_gwas_catalog_sig_filtered.annotated.tsv | tail -n +2 | grep -v 'coding_exon' | grep -Fwf prostate_noncoding_germline.tmp.list | awk '{print $0 "\tprostate"}' > prostate_noncoding_snps.tmp.tsv

#New Prostate Section
cut -f5,11,26,27 new_pca_only_gwas_catalog_sig_filtered.annotated.tsv | tail -n +2 | grep -v 'coding_exon' | grep -Fwf prostate_noncoding_germline.tmp.list | awk '{print $2 "\t" $1 "\t" $3 "\t" $4 "\tprostate"}' > new_prostate_noncoding_snps.tmp.tsv

#combine all files
cat breast_noncoding_snps.tmp.tsv colorectal_noncoding_snps.tmp.tsv lung_noncoding_snps.tmp.tsv kidney_noncoding_snps.tmp.tsv prostate_noncoding_snps.tmp.tsv new_prostate_noncoding_snps.tmp.tsv | grep -v ',' | awk 'BEGIN {FS="\t"; OFS="\t"} {split($1, a, "-"); print $3 "\t" $2 "\t" a[2] "\t" $5}' | awk 'BEGIN {FS="\t"; OFS="\t"} {split($1, a, ":"); print "chr"a[1] "\t" a[2]-1 "\t" a[2] "\t" $2 "\t" $3 "\t" $4 "\t" $5}' | sort | uniq > pancan_gwas_catalog_sig_filtered.tmp.tsv

grep 'HLA-' pancan_gwas_catalog_sig_filtered.tmp.tsv > HLA.tmp.tsv
grep -v '-' pancan_gwas_catalog_sig_filtered.tmp.tsv > pancan_gwas_catalog_sig_filtered.tmp2.tsv
cat HLA.tmp.tsv >> pancan_gwas_catalog_sig_filtered.tmp2.tsv

sort -V pancan_gwas_catalog_sig_filtered.tmp2.tsv | uniq > pancan_gwas_catalog_sig_filtered.tsv

#remove all tmp files
rm *tmp*
