#!/bin/bash
COSMIC_TABLE="/Users/noah/Desktop/DFCI_Data/gsi/data/COSMIC.CGC.10_26_23.tsv"
CELLCHAT_TABLE="/Users/noah/Desktop/DFCI_Data/gsi/data/cellchat_db_formatted.csv"
VALab_germline_somatic_table="/Users/noah/Desktop/DFCI_Data/gsi/data/VALab_germline_somatic.tsv"
output_file="/Users/noah/Desktop/DFCI_Data/gsi/data/VALab_germline_somatic.annotated.tsv"

# Function to check ONC status from COSMIC.CGC.10_26_23.tsv
get_onc_status() {
    local gene="$1"
    local oncogene="$(grep -m 1 "^$gene" "$COSMIC_TABLE" | cut -d $'\t' -f 15)"
    
    if grep -q "oncogene" <<< "$oncogene" && grep -q "TSG" <<< "$oncogene"; then
        echo "ONC-TSG"
    elif grep -q "oncogene" <<< "$oncogene"; then
        echo "ONC"
    elif grep -q "TSG" <<< "$oncogene"; then
        echo "TSG"
    else
        echo "UNK"
    fi
}

# Function to check Signaling Scenario from cellchat_db_formatted.csv
get_signaling_scenario() {
    local gene="$1"
    
    local signaling_scenario=""
    
    if grep -q -F "$gene" <(cut -d ',' -f4 "$CELLCHAT_TABLE");then
	signaling_scenario="Ligand"
    elif grep -q -F "$gene" <(cut -d ',' -f5 "$CELLCHAT_TABLE");then
	signaling_scenario="Receptor"
    else
	signaling_scenario="UNK"
    fi
    
    echo "$signaling_scenario"
}

get_mutation_type() {
    local gene="$1"
    
    local mutation_type="$(grep -m 1 "^$gene" "$COSMIC_TABLE" | cut -d $'\t' -f16)"

    if test -z "$mutation_type";then
        echo "UNK"
    else
        echo "$mutation_type"
    fi
}

echo -e "cancer\tgermline_gene\tgermline_context\tgermline_onc_status\tgermline_mutations\tgermline_signaling\tsomatic_gene\tsomatic_context\tsomatic_onc_status\tsomatic_mutations\tsomatic_signaling\tcriteria" > "$output_file"
# Main script
while IFS=$'\t' read -r cancer germline_gene germline_context somatic_gene somatic_context criteria; do
    # Get ONC status for germline gene
    germline_onc_status=$(get_onc_status "$germline_gene")
    
    # Get ONC status for somatic gene
    somatic_onc_status=$(get_onc_status "$somatic_gene")
    
    # Get Signaling Scenario
    germline_signaling_scenario=$(get_signaling_scenario "$germline_gene")
    somatic_signaling_scenario=$(get_signaling_scenario "$somatic_gene")

    #Get Mutation Type
    germline_mutation_type=$(get_mutation_type "$germline_gene")
    somatic_mutation_type=$(get_mutation_type "$somatic_gene")

    # Print the updated line
    echo -e "$cancer\t$germline_gene\t$germline_context\t$germline_onc_status\t$germline_mutation_type\t$germline_signaling_scenario\t$somatic_gene\t$somatic_context\t$somatic_onc_status\t$germline_mutation_type\t$somatic_signaling_scenario\t$criteria" >> "$output_file"
done < "$VALab_germline_somatic_table" | tail -n +2

