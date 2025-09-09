#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Summarize resource usage for all tasks from single Cromwell workflow

wid=$1
EXEC_BUCKET=$2 # Optional prefix for Cromwell execution bucket

usage() {
cat << EOF

USAGE: ./summarize_workflow_resources.sh workflow_id [execution_bucket]

       Collates data from montitoring logs from a workflow executed via Cromwell

       Can optionally provide execution bucket URI
       (default: \$WORKSPACE_BUCKET/cromwell-execution/)

EOF
}


# Parse arguments
if [ $# -lt 1 ]; then
  echo -e "\nERROR: must supply workflow ID as first positional argument"
  usage
  exit
fi

if [ -z $EXEC_BUCKET ]; then
  BUCKET="$WORKSPACE_BUCKET/cromwell-execution"
fi

WRKDIR=`mktemp -d`

# Gather list of monitoring.logs
gsutil -m ls $EXEC_BUCKET/$wid/**monitoring.log \
> $WRKDIR/log.uris.list
if [ $( cat $WRKDIR/log.uris.list | wc -l ) -lt 1 ]; then
  echo -e "Error: no files named monitoring.log found in $EXEC_BUCKET/$wid/\n"
  exit 1
fi

# Localize monitoring logs while mapping to directory structure of URIs
while read uri; do
  # Recursively create all subdirectories
  # Copy log to local
done < $WRKDIR/log.uris.list

# Extract resource information from each monitoring log

# Summarize resources by task & print to stdout

# Clean up
rm -rf $WRKDIR
