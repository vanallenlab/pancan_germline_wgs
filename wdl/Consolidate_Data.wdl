 The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Noah Fields and the Dana-Farber Cancer Institute
# Contact: Noah Fields <Noah_Fieldss@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

#The purpose of this code is to count occurences of
#somatic drivers in genes where gene abbreviation is the input
#Ex: APC, VHL, BRCA1


version 1.0

task get_individuals {
  input {
	File gcs_cohort_samples_list_dir
  }
  
  command <<<
  >>>
  
  runtime {
    docker: "ubuntu:latest"
  }
  output{
		Array[String] out = read_lines(~{gcs_cohort_samples_list_dir}) 
  }
}

task consolidate_charr {
	input {
		String individual
    	String cohort
	}
    command <<<
		gsutil cp gs://dfci-g2c-inputs/~{cohort}/charr/{individual}.txt .
		cut -f2-9 < {individual}.txt > charr_data.txt
    >>>
    runtime {
    	docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
	}
    output {
		String out = read_string("charr_data.txt")
	}
}

task write_output {
	input {
    	Array[String] individuals
        Array[String] charr_data
	}
    command <<<
		write_lines("${write_output}") {
            Array[String] combined_data = zip(individuals, charr_data)
            for (Array[String] pair in combined_data) {
                "${pair[0]}\t${pair[1]}"
            }
        }
    >>>
    output {
        File output_file = stdout()
    }
    runtime {
        docker: "ubuntu:latest"
    }
}

workflow consolidate {
	input {
    	File gcs_cohort_dir
        String cohort
	}
    call get_individuals{
    	input:
			gcs_cohort_dir = gcs_cohort_dir
    }
    scatter(i in get_individuals.out){
    	call consolidate_charr{
    		input:
        		individual = i
                cohort = cohort
    	}
    }
    call write_output {
		input:
        	individuals = get_individuals.out
			charr_data = consolidate_charr.out
	}
}
