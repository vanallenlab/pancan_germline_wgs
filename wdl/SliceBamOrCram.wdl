# WDL to slice a BAM or CRAM to a prespecified set of intervals
# Contact: Ryan_Collins@dfci.harvard.edu


version 1.0


workflow SliceBamOrCram {
  input {
    File bam_or_cram
    File bam_or_cram_idx
    File ref_fa
    File regions_bed
    String suffix = "sliced"
    Boolean cram_out = true
    String gotc_docker
  }

  Boolean is_cram = sub(basename(bam_or_cram), ".*\\.", "") == "cram"
  String sample_basename = if is_cram then basename(bam_or_cram, ".cram") else basename(bam_or_cram, ".bam")
  String out_fmt = if cram_out then "cram" else "bam"
  String outfile = sample_basename + "." + suffix + "." + out_fmt
  String outfile_idx = if cram_out then outfile + ".crai" else outfile + ".bai"

  call Slice {
    input:
      bam_or_cram = bam_or_cram,
      bam_or_cram_idx = bam_or_cram_idx,
      is_cram = is_cram,
      ref_fa = ref_fa,
      regions_bed = regions_bed,
      outfile = outfile,
      outfile_idx = outfile_idx,
      cram_out = cram_out,
      docker = gotc_docker
  }
  
  output {
    File sliced_bam_or_cram = Slice.sliced_file
    File sliced_bam_or_cram_idx = Slice.sliced_file_idx
  }
}


task Slice {
  input {
    File bam_or_cram
    File bam_or_cram_idx
    Boolean is_cram
    File ref_fa
    File regions_bed
    String outfile
    String outfile_idx
    Boolean cram_out
    Int? mem_gb
    String docker
  }
  Float disk_size = ( 2 * size([bam_or_cram, ref_fa], "GB") ) + 20.0

  command <<<
    set -eu -o pipefail

    # Build commands for slicing & indexing
    slice_cmd="samtools view -L ~{regions_bed} -o ~{outfile}"
    index_cmd="samtools index"
    if [ "~{is_cram}" == "true" ]; then
      samtools_options="$samtools_options -T ~{ref_fa}"
    fi
    if [ "~{cram_out}" == "true" ]; then
      samtools_options="$samtools_options -C"
    else
      samtools_options="$samtools_options -b"
      index_cmd="$index_cmd -b"
    fi
    slice_cmd="$slice_cmd ~{bam_or_cram}"
    index_cmd="$index_cmd ~{outfile}"

    # Slice
    echo -e "Slicing input with the following command:\n$slice_cmd"
    eval $slice_cmd

    # Index
    echo -e "Indexing output with the following command:\n$index_cmd"
    eval $index_cmd
  >>>

  runtime {
    docker: docker
    memory: select_first([mem_gb, 15]) + " GB"
    disks: "local-disk " + disk_size + " HDD"
    preemptible: 3
 }

  output {
    File sliced_file = "~{outfile}"
    File sliced_file_idx = "~{outfile_idx}"
  }
}
