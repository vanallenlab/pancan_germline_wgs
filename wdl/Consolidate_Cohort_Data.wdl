# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# The purpose of this code is to consolidate data for each G2C cohort.

version 1.0

task get_individuals {
  input{
  	String cohort
  }
  
  command <<<
	gsutil cat "gs://dfci-g2c-inputs/sample-lists/~{cohort}.samples.list" | awk 'NF > 0' > samples.list
  >>>
  
  runtime {
	docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
  }
  output{
	Array[String] out = read_lines("samples.list") 
  }
}

task unzip_ploidy_tsv {
  input{
  	String ploidy_location
  }
  
  command <<<
	gsutil -m cp ~{ploidy_location} ploidy.tsv.gz
	gunzip ploidy.tsv.gz
  >>>
  
  runtime {
	docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
  }
  output{
	File out = "ploidy.tsv"
  }
}


task get_sample_data{
	input{
		String sample
		String cohort
        File ploidy_tsv
    }
	command <<<
		charr_output=$(gsutil cat gs://dfci-g2c-inputs/~{cohort}/charr/~{sample}.txt | cut -f2-9 | tail -n 1)
		if [ -z "$charr_output" ]; then
			charr_output="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
		fi
		demographics_output=$(gsutil cat gs://dfci-g2c-inputs/~{cohort}/demographics/~{sample}.txt | cut -f2-5,7-9,13-15 | sed -n '2p')
		if [ -z "$demographics_output" ]; then
			demographics_output="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
		fi
		melt_INS_count=$(gsutil cat gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/melt_~{sample}.vcf.tsv | grep 'vcf_INS_count' | cut -f2)
		if [ -z "$melt_INS_count" ] || [ "$melt_INS_count" -eq 0 ]; then
			melt_INS_count="NA"
		fi

		gsutil -m cp gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/manta_~{sample}.vcf.tsv manta_~{cohort}.tmp
		gsutil -m cp gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/wham_~{sample}.vcf.tsv wham_~{cohort}.tmp
		gsutil -m cp gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/~{sample}.raw-counts.tsv metrics_~{cohort}.tmp

		rd_median=$(grep 'rd_q50' metrics_~{cohort}.tmp | cut -f2)
		if [ -z "$rd_median" ] || [ "$rd_median" -eq 0 ]; then
			rd_median="NA"
		fi

		rd_mean=$(grep 'rd_mean' metrics_~{cohort}.tmp | cut -f2)
		if [ -z "$rd_mean" ] || [ "$rd_mean" -eq 0 ]; then
			rd_mean="NA"
		fi

		manta_DEL_count=$(grep 'vcf_DEL_count' manta_~{cohort}.tmp | cut -f2)
		if [ -z "$manta_DEL_count" ] || [ "$manta_DEL_count" -eq 0 ]; then
			manta_DEL_count="NA"
		fi

		manta_DUP_count=$(grep 'vcf_DUP_count' manta_~{cohort}.tmp | cut -f2)
		if [ -z "$manta_DUP_count" ] || [ "$manta_DUP_count" -eq 0 ]; then
			manta_DUP_count="NA"
		fi

		manta_INS_count=$(grep 'vcf_INS_count' manta_~{cohort}.tmp | cut -f2)
		if [ -z "$manta_INS_count" ] || [ "$manta_INS_count" -eq 0 ]; then
			manta_INS_count="NA"
		fi

		manta_INV_count=$(grep 'vcf_INV_count' manta_~{cohort}.tmp | cut -f2)
		if [ -z "$manta_INV_count" ] || [ "$manta_INV_count" -eq 0 ]; then
			manta_INV_count="NA"
		fi

		manta_BND_count=$(grep 'vcf_BND_count' manta_~{cohort}.tmp | cut -f2)
		if [ -z "$manta_BND_count" ] || [ "$manta_BND_count" -eq 0 ]; then
			manta_BND_count="NA"
		fi

		wham_DEL_count=$(grep 'vcf_DEL_count' wham_~{cohort}.tmp | cut -f2)
		if [ -z "$wham_DEL_count" ] || [ "$wham_DEL_count" -eq 0 ]; then
			wham_DEL_count="NA"
		fi

		wham_DUP_count=$(grep 'vcf_DUP_count' wham_~{cohort}.tmp | cut -f2)
		if [ -z "$wham_DUP_count" ] || [ "$wham_DUP_count" -eq 0 ]; then
			wham_DUP_count="NA"
		fi

		ploidy_estimate=$(awk -F'\t' '$1 == "~{cohort}" {print}'  ~{ploidy_tsv} | awk -v sample="~{sample}" '{ if ($2==sample) print }' | cut -f3-26,28-30)
		if [ -z "$ploidy_estimate" ]; then
			ploidy_estimate="NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA"
		fi
		echo -e "~{sample}\t~{cohort}\t${demographics_output}\t${charr_output}\t${rd_median}\t${rd_mean}\t${manta_DEL_count}\t${manta_DUP_count}\t${manta_INS_count}\t${manta_INV_count}\t${manta_BND_count}\t${melt_INS_count}\t${wham_DEL_count}\t${wham_DUP_count}\t${ploidy_estimate}" > ~{sample}.txt
	>>>
	runtime {
		docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
	}
	output{
		String out = read_string("~{sample}.txt")
	}
}

task write_output {
	input {
		String cohort
		Array[String] data
	}
    command <<<
		echo -e "Sample\tCohort\tgrafpop_SNPs\tgrafpop_GD1\tgrafpop_GD2\tgrafpop_GD3\tpct_EUR\tpct_AFR\tpct_ASN\tchrX_count\tchrY_count\tancestry\tHQ_HOM\tHQ_HOM_RATE\tHQ_HET\tHQ_HET_RATE\tCHARR\tMEAN_REF_AB_HOM_ALT\tHETEROZYGOSITY_RATE\tINCONSISTENT_AB_HET_RATE\trd_median\trd_mean\tmanta_DEL_count\tmanta_DUP_count\tmanta_INS_count\tmanta_INV_count\tmanta_BND_count\tmelt_INS_count\twham_DEL_count\twham_DUP_count\tchr1_CopyNumber\tchr2_CopyNumber\tchr3_CopyNumber\tchr4_CopyNumber\tchr5_CopyNumber\tchr6_CopyNumber\tchr7_CopyNumber\tchr8_CopyNumber\tchr9_CopyNumber\tchr10_CopyNumber\tchr11_CopyNumber\tchr12_CopyNumber\tchr13_CopyNumber\tchr14_CopyNumber\tchr15_CopyNumber\tchr16_CopyNumber\tchr17_CopyNumber\tchr18_CopyNumber\tchr19_CopyNumber\tchr20_CopyNumber\tchr21_CopyNumber\tchr22_CopyNumber\tchrX_CopyNumber\tchrY_CopyNumber\tmedian_coverage\twgd_score\tnondiploid_bins" > ~{cohort}.consolidated_data.tsv
		cat ~{write_lines(data)} >> ~{cohort}.consolidated_data.tsv
		gsutil -m cp ~{cohort}.consolidated_data.tsv gs://fc-secure-826914ff-6f0b-48bd-b6b1-1002dbffd5a3/
    >>>

    runtime {
		docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    }
    output{
    }
}

workflow consolidate {
	input {
		String cohort
        String ploidy_location
	}
    call get_individuals{
		input:
			cohort = cohort
    }
    call unzip_ploidy_tsv{
    	input:
        	ploidy_location = ploidy_location
    }
    scatter(i in get_individuals.out){
		call get_sample_data{
			input:
				sample = i,
				cohort = cohort,
                ploidy_tsv = unzip_ploidy_tsv.out
    	}
    }
    call write_output {
		input:
			cohort = cohort,
			data = get_sample_data.out
	}
}
