# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Helper subroutine to preprocess site benchmarking data prior to QC plotting
# Not really intended to be called as a stand-alone workflow, 
# but is called from within PlotVcfQcMetrics.wdl


version 1.0


import "QcTasks.wdl" as QcTasks
import "Utilities.wdl" as Utils


workflow PrepSiteBenchDataToPlot {
  input {
    String dataset_name
    String dataset_prefix
    Array[String] interval_set_names

    Array[Array[File]]? common_snv_ppv_beds
    Array[Array[File]]? common_indel_ppv_beds
    Array[Array[File]]? common_sv_ppv_beds

    Array[Array[File]]? ppv_by_freqs
    Array[Array[File]]? sensitivity_by_freqs

    String bcftools_docker
    String g2c_analysis_docker
  }

  Int n_interval_sets = length(select_first([interval_set_names]))

  # Make a dummy empty BED file to help massage optional inputs
  call QcTasks.MakeEmptyBenchBed {
    input:
      docker = bcftools_docker
  }

  # Process each set of evaluation intervals in parallel
  scatter ( sbi in range(n_interval_sets) ) {

    String bi_name = interval_set_names[sbi]
    String sbi_prefix = dataset_prefix + "." + bi_name
      
    # Collapse SNV BEDs
    if (defined(common_snv_ppv_beds)) {
      Array[File] sb_common_snv_preflat = if defined(common_snv_ppv_beds) 
                                          then select_first([common_snv_ppv_beds])[sbi]
                                          else [MakeEmptyBenchBed.empty_bed]
      if ( length(sb_common_snv_preflat) > 1 ) {
        call Utils.ConcatTextFiles as CollapseSiteBenchSnvs {
          input:
            shards = sb_common_snv_preflat,
            concat_command = "zcat",
            sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
            compression_command = "bgzip -c",
            input_has_header = true,
            output_filename = sbi_prefix + ".common_snvs.bed.gz",
            docker = bcftools_docker
        }
      }
      File sb_common_snv_flat = select_first([CollapseSiteBenchSnvs.merged_file,
                                              sb_common_snv_preflat[0]])
    }
      
    # Collapse indel BEDs
    if (defined(common_indel_ppv_beds)) {
      Array[File] sb_common_indel_preflat = if defined(common_indel_ppv_beds) 
                                            then select_first([common_indel_ppv_beds])[sbi]
                                            else [MakeEmptyBenchBed.empty_bed]
      if ( length(sb_common_indel_preflat) > 1 ) {
        call Utils.ConcatTextFiles as CollapseSiteBenchIndels {
          input:
            shards = sb_common_indel_preflat,
            concat_command = "zcat",
            sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
            compression_command = "bgzip -c",
            input_has_header = true,
            output_filename = sbi_prefix + ".common_indels.bed.gz",
            docker = bcftools_docker
        }
      }
      File sb_common_indel_flat = select_first([CollapseSiteBenchIndels.merged_file,
                                                sb_common_indel_preflat[0]])
    }
      
    # Collapse SV BEDs
    if (defined(common_sv_ppv_beds)) {
      Array[File] sb_common_sv_preflat = if defined(common_sv_ppv_beds) 
                                         then select_first([common_indel_ppv_beds])[sbi]
                                         else [MakeEmptyBenchBed.empty_bed]
      if ( length(sb_common_sv_preflat) > 1 ) {
        call Utils.ConcatTextFiles as CollapseSiteBenchSvs {
          input:
            shards = sb_common_sv_preflat,
            concat_command = "zcat",
            sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
            compression_command = "bgzip -c",
            input_has_header = true,
            output_filename = sbi_prefix + ".common_svs.bed.gz",
            docker = bcftools_docker
        }
      }
      File sb_common_sv_flat = select_first([CollapseSiteBenchSvs.merged_file,
                                                sb_common_sv_preflat[0]])
    }

    # Collapse compressed PPV tsvs
    if (defined(ppv_by_freqs)) {
      call QcTasks.SumCompressedDistribs as SumSiteBenchPpvByAf {
        input:
          distrib_tsvs = select_first([ppv_by_freqs])[sbi],
          n_key_columns = 3,
          out_prefix = sbi_prefix + ".ppv_by_freq",
          g2c_analysis_docker = g2c_analysis_docker
      }
    }

    # Collapse compressed sensitivity tsvs
    if (defined(sensitivity_by_freqs)) {
      call QcTasks.SumCompressedDistribs as SumSiteBenchSensByAf {
        input:
          distrib_tsvs = select_first([sensitivity_by_freqs])[sbi],
          n_key_columns = 3,
          out_prefix = sbi_prefix + ".sensitivity_by_freq",
          g2c_analysis_docker = g2c_analysis_docker
      }
    }
  }

  # Further collapse SNVs across interval sets
  if ( defined(common_snv_ppv_beds) ) {
    Array[File] sb_common_snv_flat_array = if defined(sb_common_snv_flat)
                                           then select_all(sb_common_snv_flat)
                                           else [MakeEmptyBenchBed.empty_bed]
    if ( length(sb_common_snv_flat_array) > 1 ) {
      call Utils.ConcatTextFiles as CollapseSiteBenchSnvsUnion {
        input:
          shards = sb_common_snv_flat_array,
          concat_command = "zcat",
          sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
          compression_command = "bgzip -c",
          input_has_header = true,
          output_filename = dataset_prefix + ".merged.common_snvs.bed.gz",
          docker = bcftools_docker
      }
    }
    File sb_common_snv_flat_union = select_first([CollapseSiteBenchSnvsUnion.merged_file,
                                                  sb_common_snv_flat_array[0]])
  }

  # Further collapse indels across interval sets
  if ( defined(common_indel_ppv_beds) ) {
    Array[File] sb_common_indel_flat_array = if defined(sb_common_indel_flat)
                                             then select_all(sb_common_indel_flat)
                                             else [MakeEmptyBenchBed.empty_bed]
    if ( length(sb_common_indel_flat_array) > 1 ) {
      call Utils.ConcatTextFiles as CollapseSiteBenchIndelsUnion {
          input:
            shards = sb_common_indel_flat_array,
            concat_command = "zcat",
            sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
            compression_command = "bgzip -c",
            input_has_header = true,
            output_filename = dataset_prefix + ".merged.common_indels.bed.gz",
            docker = bcftools_docker
      }
    }
    File sb_common_indel_flat_union = select_first([CollapseSiteBenchIndelsUnion.merged_file,
                                                    sb_common_indel_flat_array[0]])
  }

  # Further collapse SVs across interval sets
  if ( defined(common_sv_ppv_beds) ) {
    Array[File] sb_common_sv_flat_array = if defined(sb_common_sv_flat)
                                          then select_all(sb_common_sv_flat)
                                          else [MakeEmptyBenchBed.empty_bed]
    if ( length(sb_common_sv_flat_array) > 1 ) {
      call Utils.ConcatTextFiles as CollapseSiteBenchSvsUnion {
          input:
            shards = sb_common_sv_flat_array,
            concat_command = "zcat",
            sort_command = "sort -Vk1,1 -k2,2n -k3,3n",
            compression_command = "bgzip -c",
            input_has_header = true,
            output_filename = dataset_prefix + ".merged.common_svs.bed.gz",
            docker = bcftools_docker
      }
    }
    File sb_common_sv_flat_union = select_first([CollapseSiteBenchSvsUnion.merged_file,
                                                 sb_common_sv_flat_array[0]])

    # For convenience, also concatenate collapsed + union BEDs into single arrays
    Array[File] sb_common_snv_all = select_all(flatten([sb_common_snv_flat, [sb_common_snv_flat_union]]))
    Array[File] sb_common_indel_all = select_all(flatten([sb_common_indel_flat, [sb_common_indel_flat_union]]))
    Array[File] sb_common_sv_all = select_all(flatten([sb_common_sv_flat, [sb_common_sv_flat_union]]))

    # For even more convenience, dump all output arrays into a .json that can 
    # be parsed by downstream plotting tasks
    call MakeOutputMenu {
      input:
        snv_ppv_bed_uris = sb_common_snv_flat,
        indel_ppv_bed_uris = sb_common_indel_flat,
        sv_ppv_bed_uris = sb_common_sv_flat,
        snv_ppv_union_uri = sb_common_snv_flat_union,
        indel_ppv_union_uri = sb_common_indel_flat_union,
        sv_ppv_union_uri = sb_common_sv_flat_union,
        snv_ppv_beds_plus_union_uris = sb_common_snv_all,
        indel_ppv_beds_plus_union_uris = sb_common_indel_all,
        sv_ppv_beds_plus_union_uris = sb_common_sv_all,
        ppv_by_af_uris = SumSiteBenchPpvByAf.merged_distrib,
        sens_by_af_uris = SumSiteBenchSensByAf.merged_distrib,
        out_json = dataset_prefix + ".site_benchmarking_plot_inputs.json"
    }
  }

  output {
    File plot_files_json = select_first([MakeOutputMenu.json])

    Array[File?] snv_ppv_beds = sb_common_snv_flat
    Array[File?] indel_ppv_beds = sb_common_indel_flat
    Array[File?] sv_ppv_beds = sb_common_sv_flat

    File? snv_ppv_union = sb_common_snv_flat_union
    File? indel_ppv_union = sb_common_indel_flat_union
    File? sv_ppv_union_union = sb_common_sv_flat_union

    Array[File]? snv_ppv_beds_plus_union = select_first([sb_common_snv_all, []])
    Array[File]? indel_ppv_beds_plus_union = select_first([sb_common_indel_all, []])
    Array[File]? sv_ppv_beds_plus_union = select_first([sb_common_sv_all, []])
  
    Array[File?] ppv_by_af = SumSiteBenchPpvByAf.merged_distrib
    Array[File?] sens_by_af = SumSiteBenchSensByAf.merged_distrib
  }
}


# Dump the URIs of all prepared site benchmarking files to a .json
task MakeOutputMenu {
  input {
    Array[String?] snv_ppv_bed_uris
    Array[String?] indel_ppv_bed_uris
    Array[String?] sv_ppv_bed_uris

    String? snv_ppv_union_uri
    String? indel_ppv_union_uri
    String? sv_ppv_union_uri

    Array[String]? snv_ppv_beds_plus_union_uris
    Array[String]? indel_ppv_beds_plus_union_uris
    Array[String]? sv_ppv_beds_plus_union_uris

    Array[String?] ppv_by_af_uris
    Array[String?] sens_by_af_uris

    String out_json
  }

  String snv_ppv_all_json = write_json(select_first([snv_ppv_beds_plus_union_uris, []]))
  String indel_ppv_all_json = write_json(select_first([indel_ppv_beds_plus_union_uris, []]))
  String sv_ppv_all_json = write_json(select_first([sv_ppv_beds_plus_union_uris, []]))
  String ppv_by_af_json = write_json(select_first([ppv_by_af_uris, []]))
  String sens_by_af_json = write_json(select_first([sens_by_af_uris, []]))

  command <<<
    set -eu -o pipefail

    cat > input_uris.json <<EOF
{
  "snv_ppv_beds": ~{write_json(snv_ppv_bed_uris)},
  "indel_ppv_beds": ~{write_json(indel_ppv_bed_uris)},
  "sv_ppv_beds": ~{write_json(sv_ppv_bed_uris)},
  "snv_ppv_union": "~{snv_ppv_union_uri}",
  "indel_ppv_union": "~{indel_ppv_union_uri}",
  "sv_ppv_union_union": "~{sv_ppv_union_uri}",
  "snv_ppv_all": ~{snv_ppv_all_json},
  "indel_ppv_all": ~{indel_ppv_all_json},
  "sv_ppv_all": ~{sv_ppv_all_json},
  "ppv_by_af": ~{ppv_by_af_json},
  "sens_by_af": ~{sens_by_af_json}
}
EOF

    python3 <<CODE

    import json

    with open("input_uris.json") as f:
      raw = json.load(f)

    data = {
      "ppv": {
        "snv": {
          "flat": raw["snv_ppv_beds"],
          "union": raw["snv_ppv_union"],
          "all": raw["snv_ppv_all"]
        },
        "indel": {
          "flat": raw["indel_ppv_beds"],
          "union": raw["indel_ppv_union"],
          "all": raw["indel_ppv_all"]
        },
        "sv": {
          "flat": raw["sv_ppv_beds"],
          "union": raw["sv_ppv_union_union"],
          "all": raw["sv_ppv_all"]
        }
      },
      "af_metrics": {
        "ppv_by_af": raw["ppv_by_af"],
        "sens_by_af": raw["sens_by_af"]
      }
    }

    with open("~{out_json}", "w") as out:
        json.dump(data, out, indent=2)   

    CODE
  >>>

  output {
    File json = "~{out_json}"
  }

  runtime {
    docker: "python:3.9-slim"
    memory: "1.7 GB"
    cpu: 1
    disks: "local-disk 20 HDD"
    preemptible: 3
    max_retries: 1
  }
}
