#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# General helper bash functions for the All of Us Researcher Workbench


# Clean up intermediate files generated by Cromwell and collected by check_status
cleanup_garbage() {
  dt_fmt=$( date '+%m_%d_%Y.%Hh%Mm%Ss' )
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

