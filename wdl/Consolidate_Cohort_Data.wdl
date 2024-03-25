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


task get_sample_data{
	input{
		String sample
		String cohort
    }
	command <<<
		charr_output=$(gsutil cat gs://dfci-g2c-inputs/~{cohort}/charr/~{sample}.txt | cut -f2-9 | tail -n 1)
		demographics_output=$(gsutil cat gs://dfci-g2c-inputs/~{cohort}/demographics/~{sample}.txt | cut -f2-5,7-9,13-15 | sed -n '2p')
		melt_INS_count=$(gsutil cat gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/melt_~{sample}.vcf.tsv | grep 'vcf_INS_count' | cut -f2)
		gsutil -m cp gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/manta_~{sample}.vcf.tsv manta_~{cohort}.tmp
		gsutil -m cp gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/wham_~{sample}.vcf.tsv wham_~{cohort}.tmp
		gsutil -m cp gs://dfci-g2c-inputs/~{cohort}/gatk-sv/metrics/~{sample}.raw-counts.tsv metrics_~{cohort}.tmp
		rd_median=$(grep 'rd_q50' metrics_~{cohort}.tmp | cut -f2)
		rd_mean=$(grep 'rd_mean' metrics_~{cohort}.tmp | cut -f2)
		manta_DEL_count=$(grep 'vcf_DEL_count' manta_~{cohort}.tmp | cut -f2)
		manta_DUP_count=$(grep 'vcf_DUP_count' manta_~{cohort}.tmp | cut -f2)
		manta_INS_count=$(grep 'vcf_INS_count' manta_~{cohort}.tmp | cut -f2)
		manta_INV_count=$(grep 'vcf_INV_count' manta_~{cohort}.tmp | cut -f2)
		manta_BND_count=$(grep 'vcf_BND_count' manta_~{cohort}.tmp | cut -f2)
		wham_DEL_count=$(grep 'vcf_DEL_count' wham_~{cohort}.tmp | cut -f2)
		wham_DUP_count=$(grep 'vcf_DUP_count' wham_~{cohort}.tmp | cut -f2)
		echo -e "~{sample}\t~{cohort}\t${demographics_output}\t${charr_output}\t${rd_median}\t${rd_mean}\t${manta_DEL_count}\t${manta_DUP_count}\t${manta_INS_count}\t${manta_INV_count}\t${manta_BND_count}\t${melt_INS_count}\t${wham_DEL_count}\t${wham_DUP_count}" > ~{sample}.txt
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
		echo -e "Sample\tCohort\tgrafpop_SNPs\tgrafpop_GD1\tgrafpop_GD2\tgrafpop_GD3\tpct_EUR\tpct_AFR\tpct_ASN\tchrX_count\tchrY_count\tancestry\tHQ_HOM\tHQ_HOM_RATE\tHQ_HET\tHQ_HET_RATE\tCHARR\tMEAN_REF_AB_HOM_ALT\tHETEROZYGOSITY_RATE\tINCONSISTENT_AB_HET_RATE\trd_median\trd_mean\tmanta_DEL_count\tmanta_DUP_count\tmanta_INS_count\tmanta_INV_count\tmanta_BND_count\tmelt_INS_count\twham_DEL_count\twham_DUP_count" > ~{cohort}.consolidated_data.tsv
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
	}
    call get_individuals{
		input:
			cohort = cohort
    }
    scatter(i in get_individuals.out){
		call get_sample_data{
			input:
				sample = i,
				cohort = cohort
    	}
    }
    call write_output {
		input:
			cohort = cohort,
			data = get_sample_data.out
	}
}
