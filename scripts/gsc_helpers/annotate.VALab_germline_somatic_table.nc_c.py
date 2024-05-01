# Read nc_c.map and create a mapping of germline genes to list of (chr:pos, risk_allele) tuples
gene_map = {}
with open("nc_c.map", "r") as map_file:
    for line in map_file:
        parts = line.strip().split("\t")
        if len(parts) == 3:
            chrom_pos, risk_allele, gene = parts
            if gene in gene_map:
                gene_map[gene].append((chrom_pos, risk_allele))
            else:
                gene_map[gene] = [(chrom_pos, risk_allele)]

# Read VALab_germline_somatic.tsv, filter rows, and create output
output = []
with open("VALab_germline_somatic.tsv", "r") as input_file:
    # Skip header
    next(input_file)
    for line in input_file:
        parts = line.strip().split("\t")
        if len(parts) >= 6 and parts[2] == "noncoding" and parts[4] == "coding":
            cancer_type, germline_gene, _, somatic_gene, _, criteria = parts
            if germline_gene in gene_map:
                chrom_pos_risk_alleles = gene_map[germline_gene]
                for chrom_pos, risk_allele in chrom_pos_risk_alleles:
                    output.append(f"{cancer_type}\t{chrom_pos}-{risk_allele}\t{germline_gene}\tnoncoding\t{somatic_gene}\tcoding\t{criteria}")

# Write output to a new file
with open("output_file.txt", "w") as output_file:
    output_file.write("cancer\tchr:pos-risk_allele\tgermline_gene\tgermline_context\tsomatic_gene\tsomatic_context\tcriteria\n")
    for line in output:
        output_file.write(line + "\n")
