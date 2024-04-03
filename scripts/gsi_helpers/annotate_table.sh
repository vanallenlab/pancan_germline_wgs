#!/bin/bash

# Function to check ONC status from COSMIC.CGC.10_26_23.tsv
get_onc_status() {
    local gene="$1"
    local oncogene="$(grep -m 1 "^$gene" COSMIC.CGC.10_26_23.tsv | cut -d $'\t' -f 15)"
    
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
    
    if grep -q -F "$gene" <(cut -d ',' -f4 cellchat_db_formatted.csv);then
	signaling_scenario="Ligand"
    elif grep -q -F "$gene" <(cut -d ',' -f5 cellchat_db_formatted.csv);then
	signaling_scenario="Receptor"
    else
	signaling_scenario="UNK"
    fi
    
    echo "$signaling_scenario"
}

echo -e "cancer\tgermline_gene\tgermline_context\tgermline_onc_status\tgermline_signaling\tsomatic_gene\tsomatic_context\tsomatic_onc_status\tsomatic_signaling\tcriteria"
# Main script
while IFS=$'\t' read -r cancer germline_gene germline_context somatic_gene somatic_context criteria; do
    # Get ONC status for germline gene
    germline_onc_status=$(get_onc_status "$germline_gene")
    
    # Get ONC status for somatic gene
    somatic_onc_status=$(get_onc_status "$somatic_gene")
    
    # Get Signaling Scenario
    germline_signaling_scenario=$(get_signaling_scenario "$germline_gene")
    somatic_signaling_scenario=$(get_signaling_scenario "$somatic_gene")
    # Print the updated line
    echo -e "$cancer\t$germline_gene\t$germline_context\t$germline_onc_status\t$germline_signaling_scenario\t$somatic_gene\t$somatic_context\t$somatic_onc_status\t$somatic_signaling_scenario\t$criteria"
done < VALab_germline_somatic.tsv | tail -n +2

