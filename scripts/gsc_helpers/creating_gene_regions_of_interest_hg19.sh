# Example code for creating gene regions of interest and subsetting AoU VCF


# The below can all be run on your local machine from the terminal

# Download MANE GTF:
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_45/GRCh37_mapping/gencode.v45lift37.annotation.gtf.gz
gsc_table="/Users/noah/Desktop/DFCI_Data/gsi/data/VALab_germline_somatic.tsv"
cut -f2,3 "$gsc_table" | grep -v 'noncoding' | cut -f1 | sort | uniq > coding_genes.list

# Subset MANE to genes of interest
zcat gencode.v45lift37.annotation.gtf.gz \
| fgrep -wf coding_genes.list \
| bgzip -c \
> gencode.v45lift37.annotation.subsetted.gtf.gz

# Extract coordinates from subsetted GTF, pad by ±2kb to be conservative, and use BEDTools to merge these coordinates into unique regions
zcat gencode.v45lift37.annotation.subsetted.gtf.gz \
| awk -v FS="\t" -v OFS="\t" -v buffer=2000 '{ print $1, $4-buffer, $5+buffer }' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
| bgzip -c \
> coding_genes.coordinates.bed.gz

rm coding_genes.list

