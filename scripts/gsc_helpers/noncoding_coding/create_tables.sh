# This file is meant to organize for logistic regression for noncoding-coding convergences
germline_somatic_file=$1
non_coding_chr=$2
non_coding_pos=$3
somatic_gene=$4
final_output_file="${non_coding_chr}.${non_coding_pos}.${somatic_gene}.tsv"
# Here is where we compile lists
grep "$somatic_gene" "$germline_somatic_file" | grep "$non_coding_chr   $non_coding_pos" | grep 'HET' | cut -f1 | sort | uniq> apc_patients_het.list

grep "$somatic_gene" "$germline_somatic_file" | grep "$non_coding_chr   $non_coding_pos" | grep 'HOM' | cut -f1 | sort | uniq> apc_patients_risk_hom.list
cut -f1 "$germline_somatic_file" | sort | uniq > patients.list

grep "$somatic_gene" "$germline_somatic_file" | cut -f1 | sort | uniq > apc_patients.list
grep -vFwf apc_patients_het.list "$germline_somatic_file" | grep -vFwf apc_patients_risk_hom.list | cut -f1 | sort | uniq > apc_patients_safe_hom.list
# Define input files
patients_file="patients.list"
apc_patients_file="apc_patients_risk_hom.list"
output_file="risk_hom_status.list"
tmp_file=$(mktemp)
# Loop through patients.list and check if each patient ID is in apc_patients_file

# Create an array of patient IDs from apc_patients_file
while read -r patient_id; do
    if grep -w "^$patient_id$" "$apc_patients_file" > "$tmp_file"; then
        echo -e "1" >> "$output_file"  # Patient is in apc_patients_file, write "1"
    else
        echo -e "0" >> "$output_file"  # Patient is not in apc_patients_file, write "0"
    fi
done < "$patients_file"

# Define input files
apc_patients_file="apc_patients_het.list"
output_file="het_status.list"

# Create an array of patient IDs from apc_patients_file
while read -r patient_id; do
    if grep -w "^$patient_id$" "$apc_patients_file" > "$tmp_file"; then
        echo -e "1" >> "$output_file"  # Patient is in apc_patients_file, write "1"
    else
        echo -e "0" >> "$output_file"  # Patient is not in apc_patients_file, write "0"
    fi
done < "$patients_file"


# Define input files
apc_patients_file="apc_patients_safe_hom.list"
output_file="safe_hom_status.list"

# Create an array of patient IDs from apc_patients_file
while read -r patient_id; do
    if grep -w "^$patient_id$" "$apc_patients_file" > "$tmp_file"; then
        echo -e "1" >> "$output_file"  # Patient is in apc_patients_file, write "1"
    else
        echo -e "0" >> "$output_file"  # Patient is not in apc_patients_file, write "0"
    fi
done < "$patients_file"

# Define input files
apc_patients_file="apc_patients.list"
output_file="driver.list"

# Create an array of patient IDs from apc_patients_file
while read -r patient_id; do
    if grep -w "^$patient_id$" "$apc_patients_file" > "$tmp_file"; then
        echo -e "1" >> "$output_file"  # Patient is in apc_patients_file, write "1"
    else
        echo -e "0" >> "$output_file"  # Patient is not in apc_patients_file, write "0"
    fi
done < "$patients_file"

echo "patients  risk_hom        het     safe_hom        driver_status" > $final_output_file
paste patients.list risk_hom_status.list het_status.list safe_hom_status.list driver.list >> $final_output_file
rm  patients.list risk_hom_status.list het_status.list safe_hom_status.list driver.list
rm apc_patients.list apc_patients_safe_hom.list apc_patients_het.list apc_patients_risk_hom.list
