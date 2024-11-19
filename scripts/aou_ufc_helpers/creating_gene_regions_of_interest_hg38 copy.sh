# Example code for creating gene regions of interest and subsetting AoU VCF

# Given the following input files: 
# 1. list (flat .txt file) of gene symbols of relevance: riaz_genes.list
# 2. list (flat .txt file) of sample IDs to retain: keep_samples.list
# 3. VCF of all variants from AoU: aou.vcf.gz (must also have .tbi index present locally)

# The below can all be run on your local machine from the terminal

# Download MANE GTF:
wget https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/release_1.3/MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz

# Subset MANE to genes of interest
zcat MANE.GRCh38.v1.3.ensembl_genomic.gtf.gz \
| fgrep -wf riaz_genes.list \
| bgzip -c \
> MANE.GRCh38.v1.3.ensembl_genomic.subsetted.gtf.gz

# Extract coordinates from subsetted GTF, pad by ±2kb to be conservative, and use BEDTools to merge these coordinates into unique regions
zcat MANE.GRCh38.v1.3.ensembl_genomic.subsetted.gtf.gz \
| awk -v FS="\t" -v OFS="\t" -v buffer=2000 '{ print $1, $4-buffer, $5+buffer }' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bgzip -c \
> riaz_genes.coordinates.bed.gz

# Copy gene coordinates to your AoU workspace bucket
# You can find this by running "echo $WORKSPACE_BUCKET" within the AoU workbench terminal
#gsutil cp riaz_genes.coordinates.bed.gz gs://my_workspace_bucket/

# Subset variants to minimal set of interest for our purposes
# Note that this must be run in the AoU workspace terminal, not on your local machine
bcftools view \
  -O z -o aou.subsetted.vcf.gz \
  --regions-file riaz_genes.coordinates.bed.gz \
  --samples-file keep_samples.list \
  --min-ac 1 \
  aou.vcf.gz

# Index output VCF for future queries
tabix -p vcf -f aou.subsetted.vcf.gz