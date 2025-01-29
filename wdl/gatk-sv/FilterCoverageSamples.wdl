# Restrict coverage matrix and median coverage file to include only a specified set of samples

version 1.0

workflow FilterCoverageSamples {
  input {
    File? bincov
    File? median_cov
    File keep_samples_list
    String sv_base_mini_docker
  }

  call FilterCovFiles {
    input:
      bincov = bincov,
      median_cov = median_cov,
      keep_samples_list = keep_samples_list,
      docker = sv_base_mini_docker
  }

  output {
    File? filtered_bincov = FilterCovFiles.filtered_bincov
    File? filtered_bincov_idx = FilterCovFiles.filtered_bincov_idx
    File? filtered_median_cov = FilterCovFiles.filtered_median_cov
  }
}

task FilterCovFiles {
  input {
    File? bincov
    File? median_cov
    File keep_samples_list
    String docker
  }

  String bincov_out = if defined(bincov) 
    then basename(select_first([bincov]), ".txt.gz") + ".subset_samples.txt.gz"
    else ""
  String bincov_idx_out = if defined(bincov)
    then basename(select_first([bincov]), ".txt.gz") + ".subset_samples.txt.gz.tbi"
    else ""
  String medcov_out = if defined(median_cov)
    then basename(select_first([median_cov]), ".txt.gz") + ".subset_samples.txt.gz"
    else ""

  Int disk_gb = ceil(3 * size(select_all([bincov, median_cov]), "GB")) + 20

  command <<<
    set -eu -o pipefail

    # Filter bincov, if provided
    if [ ~{defined(bincov)} == "true" ]; then
      zcat ~{bincov} | head -n1 | sed 's/\t/\n/g' > bincov.header || true
      drop_idxs=$( awk -v OFS="\t" '{ print $1, NR }' bincov.header \
                   | sed '1,3d' \
                   | fgrep -wvf ~{keep_samples_list} \
                   | awk -v FS="\t" '{ print $2 }' \
                   | paste -s -d, )
      zcat ~{bincov} | cut --complement -f"$drop_idxs" | bgzip -c > ~{bincov_out}
      tabix -p bed ~{bincov_out}
    fi

    # Filter medial coverage, if provided
    if [ ~{defined(median_cov)} == "true" ]; then
      head -n1 ~{median_cov} | sed 's/\t/\n/g' > medcov.header || true
      drop_idxs=$( awk -v OFS="\t" '{ print $1, NR }' medcov.header \
                   | fgrep -wvf ~{keep_samples_list} \
                   | awk -v FS="\t" '{ print $2 }' \
                   | paste -s -d, )
      cut --complement -f"$drop_idxs" ~{median_cov} > ~{medcov_out}
    fi
  >>>

  output {
    File? filtered_bincov = bincov_out
    File? filtered_bincov_idx = bincov_idx_out
    File? filtered_median_cov = medcov_out
  }

  runtime {
    cpu: 2
    memory: "3.5 GiB"
    disks: "local-disk " + disk_gb + " HDD"
    docker: docker
    preemptible: 3
    maxRetries: 1
  }
}

