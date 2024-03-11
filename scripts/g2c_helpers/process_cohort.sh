#!/bin/bash

while getopts "c:" flag; do
    case $flag in
        c) cohort=$OPTARG ;;
        *) echo "Usage: $0 -c <cohort>"; exit 1 ;;
    esac
done

if [ -z "$cohort" ]; then
    echo "Cohort name is required. Use -c flag."
    exit 1
fi

# Download the sample list file and store it as an array
samples=()
while IFS= read -r sample; do
    samples+=("$sample")
done < <(gsutil cat "gs://dfci-g2c-inputs/sample-lists/${cohort}.samples.list")
output_file="${cohort}_consolidated_data.tsv"
echo -e "Sample\tCohort\tgrafpop_SNPs\tgrafpop_GD1\tgrafpop_GD2\tgrafpop_GD3\tpct_EUR\tpct_AFR\tpct_ASN\tchrX_count\tchrY_count\tancestry\tHQ_HOM\tHQ_HOM_RATE\tHQ_HET\tHQ_HET_RATE\tCHARR\tMEAN_REF_AB_HOM_ALT\tHETEROZYGOSITY_RATE\tINCONSISTENT_AB_HET_RATE\trd_median\trd_mean\tmanta_DEL_count\tmanta_DUP_count\tmanta_INS_count\tmanta_INV_count\tmanta_BND_count\tmelt_INS_count\twham_DEL_count\twham_DUP_count" > $output_file
# Process each sample
for sample in "${samples[@]}"; do
    # Get the string output from the specified file
    charr_output=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/charr/${sample}.txt" | cut -f2-9 | tail -n 1)
    
    demographics_output=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/demographics/${sample}.txt" | cut -f2-5,7-9,13-15 | sed -n '2p')
    rd_median=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/${sample}.raw-counts.tsv" | grep 'rd_q50' | cut -f2)
    rd_mean=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/${sample}.raw-counts.tsv" | grep 'rd_mean' | cut -f2)

    manta_DEL_count=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/manta_${sample}.vcf.tsv" | grep 'vcf_DEL_count' | cut -f2)
    manta_DUP_count=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/manta_${sample}.vcf.tsv" | grep 'vcf_DUP_count' | cut -f2)
    manta_INS_count=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/manta_${sample}.vcf.tsv" | grep 'vcf_INS_count' | cut -f2)
    manta_INV_count=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/manta_${sample}.vcf.tsv" | grep 'vcf_INV_count' | cut -f2)
    manta_BND_count=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/manta_${sample}.vcf.tsv" | grep 'vcf_BND_count' | cut -f2)

    melt_INS_count=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/melt_${sample}.vcf.tsv" | grep 'vcf_INS_count' | cut -f2)
    wham_DEL_count=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/wham_${sample}.vcf.tsv" | grep 'vcf_DEL_count' |cut -f2)
    wham_DUP_count=$(gsutil cat "gs://dfci-g2c-inputs/${cohort}/gatk-sv/metrics/wham_${sample}.vcf.tsv" | grep 'vcf_DUP_count' |cut -f2)
    # Append the string output to the sample name and print to test_output.txt
    echo -e "${sample}\t${cohort}\t${demographics_output}\t${charr_output}\t${rd_median}\t${rd_mean}\t${manta_DEL_count}\t${manta_DUP_count}\t${manta_INS_count}\t${manta_INV_count}\t${manta_BND_count}\t${melt_INS_count}\t${wham_DEL_count}\t${wham_DUP_count}" >> $output_file
done

echo "Processing complete. Results saved in ${output_file}."

