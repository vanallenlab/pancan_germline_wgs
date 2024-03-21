version 1.0

import "charr.wdl" as charr

workflow aou_charr {
    input {
        File samples_file
        String vcf_directory
        String output_directory
    }

    scatter (sample_id in read_lines(samples_file)) {
        call charr.charr {
            input:
                sample_id = sample_id,
                vcf_directory = vcf_directory
        }

        call cp_file {
            input:
                input_file = charr.output_file,
                output_directory = output_directory,
                sample_id = sample_id
        }
    }
}

task cp_file {
    input {
        File input_file
        String output_directory
        String sample_id
    }

    command {
        gsutil -m cp ~{input_file} ~{output_directory}/~{sample_id}.txt
    }
    runtime{
        docker: "us.gcr.io/google.com/cloudsdktool/google-cloud-cli:latest"
    }
}
