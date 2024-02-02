# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# WDL to run VEP on one or more input VCFs


version 1.0


import "https://raw.githubusercontent.com/vanallenlab/pancan_germline_wgs/main/wdl/Utilities.wdl" as Tasks


workflow Vep {
  input {
    Array[File] vcfs
    Array[File] vcf_idxs

    File reference_fasta
    File vep_cache_tarball # VEP cache tarball downloaded from Ensembl
    Array[String?] gnomad_vcf_uris # URIs for gnomAD sites VCFs stored in Google buckets that can be remotely sliced using tabix for each shard
    Array[File?] gnomad_vcf_indexes # Indexes corresponding to gnomad_vcf_uris
    Array[String?] gnomad_infos # INFO keys to annotate from gnomad VCFs
    Array[String?] vep_remote_files # URIs for files stored in Google buckets that can be remotely sliced using tabix for each shard
    Array[File?] vep_remote_file_indexes # Indexes corresponding to vep_remote_files
    Array[File?] other_vep_files # All other files needed for VEP. These will be localized in full to each VM and moved to execution directory.

    Array[String] vep_options = [""]
    String vep_assembly = "GRCh38"
    Int vep_version = 110

    Boolean shard_vcfs = true
    Int records_per_shard = 50000
    Int max_buffered_span = 1000000
    Int remote_query_buffer = 2
    Boolean combine_output_vcfs = false
    String? cohort_prefix

    String bcftools_docker
    String vep_docker = "vanallenlab/g2c-vep:latest"
    String g2c_docker = "vanallenlab/g2c_pipeline:latest"
  }

  Int n_remote_files = length(select_all(vep_remote_files))
  Int n_gnomad_files = length(select_all(gnomad_vcf_uris))
  Boolean any_remote = n_remote_files + n_gnomad_files > 0 

  scatter ( vcf_info in zip(vcfs, vcf_idxs) ) {
    File vcf = vcf_info.left
    File vcf_idx = vcf_info.right

    if (shard_vcfs) {
      call Tasks.ShardVcf {
        input:
          vcf = vcf,
          vcf_idx = vcf_idx,
          records_per_shard = records_per_shard,
          bcftools_docker = bcftools_docker
      }

      scatter ( shard_info in zip(ShardVcf.vcf_shards, ShardVcf.vcf_shard_idxs) ) {
        call SplitRegions {
          input:
            vcf = shard_info.left,
            region_span = max_buffered_span,
            variant_buffer = remote_query_buffer,
            docker = g2c_docker
        }

        call ShardVcfByRegions {
          input:
            vcf = shard_info.left,
            vcf_idx = shard_info.right
            scatter_regions = MakeRegions.regions
            bcftools_docker = bcftools_docker
        }
      }

      Array[File] vcf_fine_shards = flatten(ShardVcfByRegions.vcf_shards)
      Array[File] vcf_fine_shard_idxs = flatten(ShardVcfByRegions.vcf_shard_idxs)
    }
    
    Array[File] vcf_shards = select_first([vcf_fine_shards, [vcf]])
    Array[File] vcf_shard_idxs = select_first([vcf_fine_shard_idxs, [vcf_idx]])

    scatter ( shard_info in zip(vcf_shards, vcf_shard_idxs) ) {

      if (any_remote) {
        call SliceRemoteFiles {
          input:
            vcf = shard_info.left,
            vcf_idx = shard_info.right,
            gnomad_vcf_uris = gnomad_vcf_uris,
            gnomad_vcf_indexes = gnomad_vcf_indexes,
            gnomad_infos = gnomad_infos,
            vep_remote_files = vep_remote_files,
            vep_remote_file_indexes = vep_remote_file_indexes,
            query_buffer = remote_query_buffer
        }
      }

      Array[File?] all_other_vep_files = flatten(select_all([SliceRemoteFiles.remote_slices, 
                                                             SliceRemoteFiles.remote_slice_idxs, 
                                                             [SliceRemoteFiles.gnomad_vcf],
                                                             [SliceRemoteFiles.gnomad_vcf_idx],
                                                             other_vep_files]))

      call RunVep {
        input:
          vcf = shard_info.left,
          vcf_idx = shard_info.right,
          reference_fasta = reference_fasta,
          vep_cache_tarball = vep_cache_tarball,
          other_vep_files = all_other_vep_files,
          gnomad_infos = select_all(gnomad_infos),
          vep_options = vep_options,
          vep_assembly = vep_assembly,
          vep_version = vep_version,
          docker = vep_docker
      }
    }

    if (shard_vcfs) {
      call Tasks.ConcatVcfs as ConcatInnerShards {
        input:
          vcfs = RunVep.annotated_vcf,
          vcf_idxs = RunVep.annotated_vcf_idx,
          out_prefix = basename(vcf, ".vcf.gz") + ".vep",
          bcftools_docker = bcftools_docker
      }
    }
    File vepped_vcf = select_first([ConcatInnerShards.merged_vcf, RunVep.annotated_vcf[0]])
    File vepped_vcf_idx = select_first([ConcatInnerShards.merged_vcf_idx, RunVep.annotated_vcf_idx[0]])
  }

  if ( combine_output_vcfs ) {
    call Tasks.ConcatVcfs as ConcatOuterShards {
      input:
        vcfs = vepped_vcf,
        vcf_idxs = vepped_vcf_idx,
        out_prefix = select_first([cohort_prefix, "all_input_vcfs"]) + ".vep.merged",
        bcftools_docker = bcftools_docker
    }
  }

  output {
    Array[File] annotated_vcfs = vepped_vcf
    Array[File] annotated_vcf_idxs = vepped_vcf_idx
    File? combined_annotated_vcf = ConcatOuterShards.merged_vcf
    File? combined_annotated_vcf_idx = ConcatOuterShards.merged_vcf_idx
  }
}


task SliceRemoteFiles {
  input {
    String vcf
    File vcf_idx
    Array[String?] gnomad_vcf_uris
    Array[File?] gnomad_vcf_indexes
    Array[String?] vep_remote_files
    Array[File?] vep_remote_file_indexes
    Array[String?] gnomad_infos

    Int query_buffer = 1000
    String gnomad_outfile_prefix = "gnomad"

    Int disk_gb = 550
    Float mem_gb = 15.5
    Int n_cpu = 4

    String docker = "us.gcr.io/broad-dsde-methods/gatk-sv/sv-base:2023-07-28-v0.28.1-beta-e70dfbd7"
  }

  Int n_remote_files = length(select_all(vep_remote_files))
  Int n_gnomad_files = length(select_all(gnomad_vcf_uris))

  command <<<
    set -eu -o pipefail

    export GCS_OAUTH_TOKEN=`gcloud auth application-default print-access-token`

    echo -e "\nSLICING QUERY REGIONS FROM VCF, REMOTELY\n"
    mv ~{vcf_idx} ./
    bcftools query --format '%CHROM\t%POS\n' ~{vcf} \
    | awk -v OFS="\t" -v buf=~{query_buffer} \
      '{ print $1, $2-buf, $2+buf }' \
    | awk -v OFS="\t" '{ if ($2<1) $2=1; print $1, $2, $3 }' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | bedtools merge -i - \
    | bgzip -c \
    > query.bed.gz
    n_queries=$( zcat query.bed.gz | wc -l )

    if [ ~{n_remote_files} -gt 0 ]; then
      echo -e "\nSLICING REMOTE ANNOTATION FILES:\n"
      mkdir remote_slices
      mv ~{sep=" " select_all(vep_remote_file_indexes)} ./

      while read uri; do
        local_name=$( basename $uri )
        echo -e "$local_name"
        if [ $n_queries -gt 0 ]; then
          tabix -h -R query.bed.gz $uri \
          | bgzip -c > remote_slices/$local_name
        else
          tabix --only-header $uri \
          | bgzip -c > remote_slices/$local_name
        fi
        tabix -s 1 -b 2 -e 2 -f remote_slices/$local_name
      done < ~{write_lines(select_all(vep_remote_files))}
    fi

    if [ ~{n_gnomad_files} -gt 0 ]; then
      echo -e "\nSLICING GNOMAD VCFs:\n"
      mkdir gnomad_slices/
      mv ~{sep=" " select_all(gnomad_vcf_indexes)} ./

      query_fmt="%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER"
      query_keys=$( cat ~{write_lines(select_all(gnomad_infos))} | awk '{ print $1"=%INFO/"$1 }' | paste -s -d\; )
      query_fmt="$query_fmt\t$query_keys\n"
      echo "Interpreted gnomAD query format as $query_fmt"

      if [ $n_queries -gt 0 ]; then
        while read uri; do
          local_name=$( basename $uri )
          echo -e "$local_name"
          tabix -H $uri | cut -f1-8 > "$local_name".header
          bcftools query \
            -R query.bed.gz \
            -f "$query_fmt" \
            $uri \
          | cat "$local_name".header - | bgzip -c \
          > gnomad_slices/"$local_name"
        done < ~{write_lines(select_all(gnomad_vcf_uris))}

        echo -e "\nMERGING GNOMAD VCFs:\n"
        bcftools concat --naive \
          -O z -o ~{gnomad_outfile_prefix}.vcf.gz \
          gnomad_slices/*.vcf.*z
      else
        tabix -H ~{select_all(gnomad_vcf_uris)[0]} \
        | bgzip -c > ~{gnomad_outfile_prefix}.vcf.gz
      fi
      tabix -p vcf -f ~{gnomad_outfile_prefix}.vcf.gz
    fi

  >>>

  output {
    Array[File?] remote_slices = glob("remote_slices/*gz")
    Array[File?] remote_slice_idxs = glob("remote_slices/*gz.tbi")
    File? gnomad_vcf = "~{gnomad_outfile_prefix}.vcf.gz"
    File? gnomad_vcf_idx = "~{gnomad_outfile_prefix}.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    bootDiskSizeGb: 30
  }
}


task RunVep {
  input {
    File vcf
    File vcf_idx

    File reference_fasta
    File vep_cache_tarball
    Array[File?] other_vep_files
    String vep_assembly

    Array[String] vep_options
    Array[String] gnomad_infos = []
    Int vep_max_sv_size = 50
    Int vep_version = 110

    Float mem_gb = 7.5
    Int n_cpu = 4
    Int? disk_gb

    String docker
  }

  String out_filename = basename(vcf, ".vcf.gz") + ".vep.vcf.gz"
  Int default_disk_gb = ceil(10 * size([vcf, vep_cache_tarball, reference_fasta], "GB")) + 50

  command <<<
    set -eu -o pipefail

    # Unpack contents of cache into $VEP_CACHE/
    # Note that $VEP_CACHE is a default ENV variable set in VEP docker
    tar -xzvf ~{vep_cache_tarball} -C $VEP_CACHE/

    # Relocate other_vep_files to execution directory
    if [ ~{defined(other_vep_files)} == "true" ]; then
      while read file; do
        mv $file ./
      done < ~{write_lines(select_all(other_vep_files))}
    fi

    # Build gnomad annotation command based on gnomad_infos
    gnomad_option=""
    if [ ~{length(gnomad_infos)} -gt 0 ]; then
      gnomad_option="--custom gnomad.vcf.gz,gnomAD,vcf,exact,0,~{sep=',' gnomad_infos}"
    fi

    vep \
      --input_file ~{vcf} \
      --format vcf \
      --output_file ~{out_filename} \
      --vcf \
      --verbose \
      --compress_output bgzip \
      --force_overwrite \
      --species homo_sapiens \
      --assembly ~{vep_assembly} \
      --max_sv_size ~{vep_max_sv_size} \
      --offline \
      --cache \
      --dir_cache $VEP_CACHE/ \
      --cache_version ~{vep_version} \
      --dir_plugins $VEP_PLUGINS/ \
      --fasta ~{reference_fasta} \
      --minimal \
      --nearest gene \
      --distance 10000 \
      --numbers \
      --hgvs \
      --no_escape \
      --symbol \
      --canonical \
      --domains \
      ~{sep=" " vep_options} \
      $gnomad_option

    tabix -f ~{out_filename}

  >>>

  output {
    File annotated_vcf = "~{out_filename}"
    File annotated_vcf_idx = "~{out_filename}.tbi"
  }

  runtime {
    docker: docker
    memory: mem_gb + " GB"
    cpu: n_cpu
    disks: "local-disk " + select_first([disk_gb, default_disk_gb]) + " HDD"
    bootDiskSizeGb: 25
    preemptible: 3
  }
}


task SplitRegions {
  input {
    File vcf
    Int region_span = 1000000
    Int variant_buffer = 100
    String docker
  }

  # Tentative
  String out_prefix = basename(basename(basename(vcf, ".gz"), ".vcf"), ".bcf")
  Int disk_gb = ceil(1.3 * size(vcf, "GB"))

  command <<<
    bcftools view -O v -o ~{out_prefix}.vcf --threads 2 ~{vcf}

    pip install --no-input numpy
  
    python3 - <<'__script__'
    import math
    from collections import deque
    import numpy as np

    in_vcf = "~{out_prefix}.vcf"
    max_span = ~{region_span}
    buffer = ~{variant_buffer}

    xrms = np.loadtxt(in_vcf, usecols=(0), dtype='U')
    pos = np.loadtxt(in_vcf, usecols=(1), dtype=np.int32)

    # unfortunately numpy sorts before returning the unique elements,
    #  so we have to unsort (i.e. re-sort the indices into order)
    Xorder = list(xrms[sorted(np.unique(xrms, return_index=True)[1])])

    L, U = np.clip(pos - buffer, 1, None), pos + buffer             # Lower, Upper buffers around each variant
    I = np.nonzero((L[1:] > U[:-1]) | (xrms[1:] != xrms[:-1]))[0]   # where buffer boundaries _dont_ overlap
    iA, iB = np.concatenate([[0], I+1]), np.concatenate([I, [-1]])  # lower, upper merged buffer indices
    X, A, B = xrms[iA], L[iA], U[iB]                                # merged buffers (i.e. subregions)
    S = (B-A).sum()                                                 # total subregion span
    N = math.ceil(S/max_span)                                       # number of regions
    T = math.ceil(S/N)                                              # target region length (for equal distribution)

    print(S, 'sized buffer to', T, 'length fragments')

    queue = deque(zip(X, A, B))
    queue.reverse() # deques pop from the end, so reverse it
    s, i, N = 0, 0, True
    with open(f"~{out_prefix}.scatter_regions.txt", 'w') as out:
      while len(queue) > 0:
        x, a, b = queue.pop()   # deque supports O(1) resizing
        if N: xl, al = x, a     # if the region is new we need to start a new record
        elif x != xl:           # if the region didn't fill and the chromosome changes...  
          out.write(f'{xl}:{al}-{bl}\t{i}\n') # we need to record what we had...
          xl, al = x, a                       # and start a new record
        if b-a > (r := T-s):    # split a subregion if it overflows a region, r is what will fit
          m = a + r             # mark the new end point to fill the region; half open intervals so m is end and next start
          queue.append((x,m,b)) # the chromosome won't change next pass, but the region will. ...
          b = m                 #   while it's unlikely for this to not happen, it's not impossible
        s += b-a                # contribute subregion measure, half open intervals -> no +1
        if (N := s >= T):       # N signals a new region/this region filled
          out.write(f'{x}:{al}-{b}\t{i}\n')   # x == xl here, so al and b are colinear
          i += 1                # increment region count
          s = 0                 # reset region measure
        bl = b                  # remember where we left off in case the chromosome changes
      if not N: out.write(f'{x}:{al}-{b}\t{i}\n') # If the region didnt fill and there's no next subregion, we need to record what we have

    __script__


  >>>

  output {
    File regions = "~{out_prefix}.scatter_regions.txt"
  }

  runtime {
    cpu: 2
    memory: "6 GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: 3
    docker: docker
  }
}


task ShardVcfByRegion {
  input {
    File vcf
    File vcf_idx
    File scatter_regions
    String bcftools_docker
  }

  String out_prefix = basename(vcf, ".vcf.gz") + ".sharded"
  Int disk_gb = ceil(3 * size(vcf, "GB"))

  command <<<
    set -eu -o pipefail

    # Make an empty shard in case the input VCF is totally empty
    bcftools view -h ~{vcf} | bgzip -c > "~{out_prefix}.0.vcf.gz"

    bcftools +scatter \
      -O z3 -o . -p "~{out_prefix}". \
      -S ~{scatter_regions} \
      ~{vcf}

    # Print all VCFs to stdout for logging purposes
    find ./ -name "*.vcf.gz"

    # Index all shards
    find ./ -name "~{out_prefix}.*.vcf.gz" \
    | xargs -I {} tabix -p vcf -f {}
  >>>

  output {
    Array[File] vcf_shards = glob("~{out_prefix}.*.vcf.gz")
    Array[File] vcf_shard_idxs = glob("~{out_prefix}.*.vcf.gz.tbi")
  }

  runtime {
    cpu: 2
    memory: "3.75 GiB"
    disks: "local-disk " + disk_gb + 20 + " HDD"
    bootDiskSizeGb: 10
    docker: bcftools_docker
    preemptible: 3
    maxRetries: 1
  }
}