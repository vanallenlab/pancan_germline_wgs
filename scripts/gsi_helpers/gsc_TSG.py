# Load COSMIC gene names into a set for efficient lookup
with open('COSMIC_TSGs.list', 'r') as cosmic_file:
    cosmic_genes = {line.strip() for line in cosmic_file}

# Define input and output file names
input_file = 'gsc_filtered.tsv'
output_file_filtered = 'gsc_filtered_cosmic.tsv'
output_file_not_in_cosmic = 'gsc_filtered_not_in_cosmic.tsv'

# Open input and output files
with open(input_file, 'r') as infile, \
        open(output_file_filtered, 'w') as outfile_filtered, \
        open(output_file_not_in_cosmic, 'w') as outfile_not_in_cosmic:

    # Write header to both output files
    header = next(infile)
    outfile_filtered.write(header)
    outfile_not_in_cosmic.write(header)

    # Iterate through each line in the input file
    for line in infile:
        parts = line.strip().split('\t')
        if len(parts) >= 6:
            gene_5th_column = parts[4]
            gene_6th_column = parts[5]

            # Check if both genes are in COSMIC
            if gene_5th_column in cosmic_genes and gene_6th_column in cosmic_genes:
                outfile_filtered.write(line)
            else:
                outfile_not_in_cosmic.write(line)

print("Filtered file with COSMIC genes:", output_file_filtered)
print("File with rows not in COSMIC:", output_file_not_in_cosmic)

