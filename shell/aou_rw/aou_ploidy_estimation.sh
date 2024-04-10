#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Shell code to run GATK-SV ploidy estimation & dosage bias scoring (module 02) on all AoU samples selected for G2C phase 1

# Note that this code is designed to be run inside the AoU Researcher Workbench
# See aou_bash_utils.sh for custom function definitions used below


#########
# SETUP #
#########

# Set up local environment
export GPROJECT="vanallen-pancan-germline-wgs"
export MAIN_WORKSPACE_BUCKET=gs://fc-secure-d21aa6b0-1d19-42dc-93e3-42de3578da45

# Prep working directory structure
for dir in cromshell cromshell/inputs cromshell/job_ids cromshell/progress \
           data data/batch_lists; do
  if ! [ -e $dir ]; then mkdir $dir; fi
done

# Copy necessary code to local disk
gsutil -m cp -r $MAIN_WORKSPACE_BUCKET/code ./
for suffix in py R; do
  find code/ -name "*.${suffix}" | xargs -I {} chmod a+x {}
done

# Source .bashrc and bash utility functions
. ~/code/refs/dotfiles/aou.rw.bashrc
. code/refs/aou_bash_utils.sh

# Install necessary packages
. code/refs/install_packages.sh

# Prep GATK-SV dependencies
cd code/wdl/gatk-sv && \
zip gatksv.dependencies.zip *.wdl && \
mv gatksv.dependencies.zip ~/ && \
cd ~

# Format local copy of Cromwell options .json to reference this workspace's storage bucket
~/code/scripts/envsubst.py \
  -i code/refs/json/aou.cromwell_options.default.json \
  -o code/refs/json/aou.cromwell_options.default.json2 && \
mv code/refs/json/aou.cromwell_options.default.json2 \
   code/refs/json/aou.cromwell_options.default.json

# Localize list of all G2C AoU samples
gsutil -m -u $GPROJECT cp \
  $MAIN_WORKSPACE_BUCKET/data/sample_info/other/all_g2c.phase1.samples.list \
  data/


#####################
# PLOIDY ESTIMATION #
#####################

# Randomly divide AoU samples into batches of ~500
code/scripts/evenSplitter.R \
  -L 500 \
  --shuffle \
  data/all_g2c.phase1.samples.list \
  data/batch_lists/AoU.ploidy_estimation.batch_

# Make list of all batches
find data/batch_lists/ -name "AoU.ploidy_estimation.batch_*" \
| xargs -I {} basename {} | sort -V \
> data/ploidy_batches.list

# Submit one workflow per batch while also checking the status of prior submissions
touch cromshell/progress/aou.ploidy_estimation.batch_status.tsv
while read bid; do

  # Check status of most recent submission, if any
  if [ -e cromshell/job_ids/$bid.ploidy_estimation.job_ids.list ]; then
    last_uuid=$( tail -n1 cromshell/job_ids/$bid.ploidy_estimation.job_ids.list )
    status=$( cromshell --no_turtle -t 120 -mc status $last_uuid \
              | jq .status | tr -d '"' | tr '[A-Z]' '[a-z]' | sed 's/[ ]+/_/g' )
  else
    status="not_started"
  fi
  echo -e "Status of batch $bid: $status"

  # Update status tracker
  fgrep -wv $bid cromshell/progress/aou.ploidy_estimation.batch_status.tsv \
  | cat - <( echo -e "$bid\t$status" ) | sort -Vk1,1 \
  >> cromshell/progress/aou.ploidy_estimation.batch_status.tsv2
  mv cromshell/progress/aou.ploidy_estimation.batch_status.tsv2 \
     cromshell/progress/aou.ploidy_estimation.batch_status.tsv

  # Only submit if status is not_started or failed
  if [ $status == "not_started" ] || [ $status == "failed" ]; then

    echo -e "Submitting ploidy estimation workflow for batch $bid"

    # Format input .json
    export BATCH=$bid
    export SAMPLES=$( awk -v ORS=", " '{ print "\""$1"\"" }' data/batch_lists/$bid | sed 's/, $//' )
    export COUNTS=$( awk -v ORS=", " -v uri_base="$MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs/aou/gatk-sv/coverage/" '{ print "\""uri_base$1".counts.tsv.gz\"" }' data/batch_lists/$bid | sed 's/, $//' )
    ~/code/scripts/envsubst.py \
      -i code/refs/json/aou.gatk_sv.02_evidence_qc.inputs.template.json \
    | sed 's/\t//g' | paste -s -d\ \
    > cromshell/inputs/$bid.ploidy_estimation.inputs.json

    # Submit job and add job ID to list of jobs for this workflow
    cmd="cromshell --no_turtle -t 120 -mc submit"
    cmd="$cmd --options-json code/refs/json/aou.cromwell_options.default.json"
    cmd="$cmd --dependencies-zip gatksv.dependencies.zip"
    cmd="$cmd code/wdl/gatk-sv/EvidenceQC.wdl"
    cmd="$cmd cromshell/inputs/$bid.ploidy_estimation.inputs.json"
    eval $cmd | jq .id | tr -d '"' \
    >> cromshell/job_ids/$bid.ploidy_estimation.job_ids.list

  fi

done < data/ploidy_batches.list






