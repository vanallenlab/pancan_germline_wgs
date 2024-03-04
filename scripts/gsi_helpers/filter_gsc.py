def filter_chrom_location(input_file, cutoff):
    # Dictionary to store counts of chromosome and location pairs
    chrom_location_counts = {}

    # Read the input TSV file
    with open(input_file, 'r') as file:
        lines = file.readlines()

    patient_cancer_dict = dict()
    num_patients = len(patient_cancer_dict.keys())
    # Count occurrences of each chromosome and location pair
    for line in lines:
        # Split the line into columns
        columns = line.strip().split('\t')

        # Extract chromosome and location
        chromosome = columns[2]
        location = columns[3]

        #Add Patient/Cancer to Dict
        patient_cancer_dict[columns[0]] = columns[1]

        # Combine chromosome and location as a tuple
        pair = (chromosome, location)
        # Increment the count for this pair in the dictionary
        chrom_location_counts[pair] = chrom_location_counts.get(pair, 0) + 1


    # Get the number of unique patients (unique values in the first column)
    num_patients = len(patient_cancer_dict.keys())

    # Filter the lines based on the counts and the cutoff
    filtered_lines = []
    filtered_test_genes_lines = []
    for line in lines:
        # Split the line into columns
        columns = line.strip().split('\t')

        # Extract chromosome and location
        chromosome = columns[2]
        location = columns[3]

        # Combine chromosome and location as a tuple
        pair = (chromosome, location)
        # Calculate the count for this pair
        pair_count = chrom_location_counts[pair]
        # Calculate the threshold count
        threshold_count = cutoff * num_patients
        # If the count is greater than the threshold, skip this line. Taking into account when chromosome == 0 for test gene 'NO_GENE'
        if pair_count > threshold_count and int(chromosome) != 0:
            continue
        # Otherwise, keep this line
        filtered_lines.append(line)
        germ_gene = columns[4]
        som_gene = columns[5].strip()
        if (germ_gene == som_gene) and (germ_gene == "APC" or germ_gene == "VHL" or germ_gene == "BRCA1"):
            filtered_test_genes_lines.append(line)

    filtered_patient_set = set(line.split('\t')[0] for line in filtered_lines)

    #We lost patients w/ Variant Filtering <5% and need to add them back for total count purposes
    count = 0 
    for patient in patient_cancer_dict.keys():
        if patient not in filtered_patient_set:
            cancer_type = patient_cancer_dict[patient]
            filtered_lines.append(f"{patient}\t{cancer_type}\t0\t0\tNO_GENE\tNO_GENE\n")

    print(len(filtered_lines))
    # Write the filtered lines to a new file
    with open('/Users/noah/Desktop/DFCI_Data/gsi/data/gsc_filtered.tsv', 'w') as outfile:
        outfile.writelines(filtered_lines)

    #with open('/Users/noah/Desktop/DFCI_Data/gsi/data/gsc_filtered_test_genes.tsv', 'w') as outfile:
        #outfile.writelines(filtered_test_genes_lines)

# Example usage:
filter_chrom_location("/Users/noah/Desktop/DFCI_Data/gsi/data/gsc.tsv", 0.05)  # Cutoff threshold of 5%
