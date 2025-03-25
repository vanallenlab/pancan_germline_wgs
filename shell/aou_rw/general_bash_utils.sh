#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# General helper bash functions for the All of Us Researcher Workbench


# Clean up intermediate files generated by Cromwell and collected by check_status
cleanup_garbage() {
  dt_fmt=$( TZ=America/New_York date '+%m_%d_%Y.%Hh%Mm%Ss' )
  garbage_uri="$WORKSPACE_BUCKET/dumpster/dfci_g2c.aou_rw.$dt_fmt.garbage"
  if [ -e ~/uris_to_delete.list ]; then
    gsutil -m cp ~/uris_to_delete.list $garbage_uri
    rm ~/uris_to_delete.list
    cat << EOF | sed -e 's/\ //g' | paste -s -d\ > ~/cromshell/inputs/empty_dumpster.$dt_fmt.inputs.json
{
  "DeleteGcpObjects.uri_list" : "$garbage_uri",
  "DeleteGcpObjects.n_cpu" : 2,
  "DeleteGcpObjects.uris_per_shard" : 100000
}
EOF
    cromshell --no_turtle -t 120 -mc submit \
      --options-json ~/code/refs/json/aou.cromwell_options.default.json \
      ~/code/wdl/pancan_germline_wgs/DeleteGcpObjects.wdl \
      ~/cromshell/inputs/empty_dumpster.$dt_fmt.inputs.json \
    | jq .id | tr -d '"' \
    >> ~/cromshell/job_ids/empty_dumpster.job_ids.list
  fi
}


# Get G2C workspace processing number from within an AoU RW terminal
get_workspace_number() {
  case "$WORKSPACE_BUCKET" in
    "gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45")
      echo 1
      ;;
    "gs://fc-secure-29075a92-7950-4778-aa20-874a75cd37bf")
      echo 2
      ;;
    "gs://fc-secure-db6876b5-357d-4339-8f79-676f87b4b2d4")
      echo 3
      ;;
    "gs://fc-secure-7b3256c8-9e62-4bb9-828b-a5d0e4b6ac31")
      echo 4
      ;;
    "gs://fc-secure-c61bace1-e344-452d-8902-49ddad15e5e7")
      echo 5
      ;;
    *)
      echo "UNKNOWN"
      ;;
  esac
}


# Concatenate a single-column text file into a Cromwell-parseable list of quoted strings
collapse_txt() {
   # Check inputs
  if [ $# -ne 1 ]; then
    echo "Must pass a single text file as input"
    return
  fi
  awk -v ORS='",' '{ print "\""$1 }' $1 | sed 's/,$//g' | awk '{ print "["$1"]" }'
}


# Find all Cromwell return code files with non-zero exit status within a given bucket prefix
# This is sometimes a useful alternative when cromshell continues to timeout
# due to a Cromwell server being overloaded
check_cromwell_return_codes() {
   # Check inputs
  if [ $# -ne 1 ]; then
    echo "Must pass a gs:// bucket prefix as only input"
    return
  fi
  bucket_prefix=$( echo $1 | sed 's/\/$//g' )
  while read uri; do
    rc=$( gsutil cat $uri )
    if [ $rc != "0" ]; then
      echo -e "$rc\t$uri"
    fi
  done < <( gsutil -m ls $bucket_prefix/**rc 2>/dev/null \
            | fgrep -v memory_retry_rc | fgrep -v cacheCopy )
  echo -e "Finished checking all return codes"
}


# Simple routine to monitor a single Cromwell workflow
# Required first positional argument is workflow ID
# Optional second positional argument is gate window (in minutes)
monitor_workflow() {
  # Check inputs
  if [ $# -lt 1 ]; then
    echo "Must provide workflow ID as first positional argument"
    return 2
  fi
  if [ $# -ge 2 ]; then
    monitor_gate=$2
  else
    monitor_gate=5
  fi

  # Endless loop
  while true; do
    echo -e "\n\n\n\n"
    date
    cromshell -t 120 --no_turtle counts -x $1
    echo -e "Waiting $monitor_gate minutes before checking again...\n"
    sleep ${monitor_gate}m
  done
}


# Update a contig-specific overrides .json for manage_chromshards.py with
# arrays of VCFs and VCF indexes staged in a google bucket
add_contig_vcfs_to_chromshard_overrides_json() {
  # Check inputs
  if [ $# -lt 4 ]; then
    cat << EOF
Error. Must provide at least four positional arguments:
  1. path to overrides.json to be updated (will be created if it doesn't exist)
  2. base URI to search for .json mapping VCF arrays. Structure of this URI
     is expected to follow the output from manage_chromshards.py; do not use 
     this function if you have a custom URI directory structure.
  3. VCF array key in .json
  4. VCF index array key in .json
  5. (Optional) path to contig list; defaults to workspace-specific contig list
EOF
    return 2
  else
    in_json=$1
    base_uri=$( echo $2 | sed 's/\/$//g' )
    vcf_key=$3
    idx_key=$4
  fi
  if [ $# -ge 5 ]; then
    contig_list=$5
  else
    contig_list=/home/jupyter/contig_lists/dfci-g2c.v1.contigs.w$WN.list
  fi

  # Make temporary directory for staging
  WRKDIR=`mktemp -d`

  # Create overrides.json if it doesn't already exist; otherwise, make a copy for updating
  if ! [ -e $in_json ]; then
    echo "{}" > $WRKDIR/updates.json
  else
    cp $in_json $WRKDIR/updates.json

  # Add each chromosome's overrides
  while read contig; do

    # Attempt to locate this chromosome's .json map of outputs
    gsutil ls $base_uri/$contig/*json > $WRKDIR/$contig.fmap.json.list
    if [ $( cat $WRKDIR/$contig.fmap.json.list | wc -l ) -lt 1 ]; then
      echo "Error: no .json was found for $contig in $base_uri/$contig/"
      rm -rf $WRKDIR && unset WRKDIR
      exit 2
    elif [ $( cat $WRKDIR/$contig.fmap.json.list | wc -l ) -gt 1 ]; then
      echo "Error: multiple .jsons were found for $contig in $base_uri/$contig/"
      rm -rf $WRKDIR && unset WRKDIR
      exit 2
    else
      contig_json=$( head -n1 $WRKDIR/$contig.fmap.json.list )
    fi

    # Write .json snippet for variable overrides for this contig
    cat << EOF > $WRKDIR/$contig.overrides.json
{
  "$contig" : {
      "CONTIG_VCFS": $( gsutil cat $contig_json | jq .$vcf_key ),
      "CONTIG_VCF_IDXS": $( gsutil cat $contig_json | jq .$idx_key )
    }
}
EOF
  
    # Update main .json
    code/scripts/update_json.py \
      -i $WRKDIR/updates.json \
      -u $WRKDIR/$contig.overrides.json \
      -o $WRKDIR/updates.json
  done < $contig_list

  # Once updates are all complete, copy updated .json back to desired location
  mv $WRKDIR/updates.json $in_json

  # Clean up
  rm -rf $WRKDIR
  unset WRKDIR
}

