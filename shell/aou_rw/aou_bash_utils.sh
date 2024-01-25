#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Helper bash functions for All of Us Researcher Workbench

# Submit workflows for one cancer type and a single workflow (gatk-hc, gatk-sv, or gcnv-pp)
submit_workflows() {
  # Check inputs
  if [ $# -ne 2 ]; then
    echo "Must specify (gatk-hc|gatk-sv|gvcf-pp) and cancer type as first two positional arguments"
  else
    wflow=$1
    cancer=$2
  fi

  # Set workflow-specific parameters
  wflow_lower=$( echo $wflow | sed 's/-/_/g' )
  wflow_nospace=$( echo $wflow_lower | sed 's/_//g' )
  case $wflow in
    "gatk-hc")
      wdl=~/code/wdl/gatk-hc/haplotypecaller-gvcf-gatk4.wdl
      inputs_json_prefix=gatk_hc
      gate_width=100
      gate_timeout=30m
      workflow_name="GATK-HC"
      sid_cram_list=~/data/cram_paths/$cancer.cram_paths.tsv
      ;;
    "gatk-sv")
      wdl=~/code/wdl/gatk-sv/GatherSampleEvidence.wdl
      inputs_json_prefix=gatk_sv_module_01
      gate_width=25
      gate_timeout=60m
      workflow_name="GATK-SV"
      sid_cram_list=~/data/cram_paths/$cancer.cram_paths.tsv
      ;;
    "gvcf-pp")
      wdl=~/code/wdl/pancan_germline_wgs/PostprocessGvcf.wdl
      inputs_json_prefix=gvcf_pp
      gate_width=100
      gate_timeout=30m
      workflow_name="gVCF postprocessing"
      awk '{ if ($2=="staged") print $1 }' \
        cromshell/progress/$cancer.gatk_hc.sample_progress.tsv \
      | fgrep -wf - data/cram_paths/$cancer.cram_paths.tsv \
      > $cancer.$wflow.sids_to_submit.list
      sid_cram_list=$cancer.$wflow.sids_to_submit.list
      ;;
  esac

  # Only submit tasks for eligible samples (depending on last known status)
  status_tsv=~/cromshell/progress/$cancer.$wflow_lower.sample_progress.tsv
  if ! [ -e $status_tsv ]; then
    echo -e "Status tracker $status_tsv not found for $cancer; skipping $workflow_name submissions for this cancer type"
    return 0
  fi
  k=0; j=0; s=0; g=0
  n=$( cat ~/data/cram_paths/$cancer.cram_paths.tsv | wc -l )
  while read sid CRAM CRAI; do
    ((k++)); ((j++))
    # Check if sample has not been launched or has failed/aborted
    status=$( awk -v FS="\t" -v sid=$sid '{ if ($1==sid) print $2 }' $status_tsv )
    if [ -z $status ] || [ $status == "not_started" ] || \
       [ $status == "failed" ] || [ $status == "aborted" ] || \
       [ $status == "unknown" ]; then

      # Format sample-specific input .json
      ((s++)); ((g++))
      export sid=$sid
      export CRAM=$CRAM
      export CRAI=$CRAI
      ~/code/scripts/envsubst.py \
        -i code/refs/json/aou.$inputs_json_prefix.inputs.template.json \
      | sed 's/\t//g' | paste -s -d\ \
      > cromshell/inputs/$sid.$wflow_nospace.inputs.json

      # Submit job and add job ID to list of jobs for this sample
      cromshell --no_turtle -t 120 -mc submit \
        --options-json code/refs/json/aou.cromwell_options.default.json \
        $wdl \
        cromshell/inputs/$sid.$wflow_nospace.inputs.json \
      | jq .id | tr -d '"' \
      >> cromshell/job_ids/$sid.$wflow_lower.job_ids.list
    fi

    # Only report on progress every 50 samples to limit terminal verbosity
    # This is one of the possible causes of large data egress warnings
    if [ $j -eq 50 ]; then
      echo -e "$k of $n $cancer samples evaluated; $s GATK-SV jobs submitted"
      j=0
    fi

    # Only submit $gate_width workflows every $gate_timeout at most
    # This is in the hopes of limiting data egress surges to not trigger security warnings
    if [ $g -eq $gate_width ]; then
      echo -e "$g $workflow_name workflows submitted recently; waiting $gate_timeout before resuming submissions..."
      g=0
      sleep $gate_timeout
    fi
  done < $sid_cram_list
  if [ $wflow == "gvcf-pp" ]; then
    rm $cancer.$wflow.sids_to_submit.list
  fi
}


# Check sample status for one cancer type and a single workflow (gatk-hc, gatk-sv, or gcnv-pp)
check_status() {
  if [ $# -ne 2 ]; then
    echo "Must specify (gatk-hc|gatk-sv|gvcf-pp) and cancer type as first two positional arguments"
  else
    wflow=$1
    cancer=$2
  fi
  status_tsv=~/cromshell/progress/$cancer.$( echo $wflow | sed 's/-/_/g' ).sample_progress.tsv
  if [ -e $status_tsv ]; then
    n=$( cat ~/sample_lists/$cancer.samples.list | wc -l )
    k=0
    j=0
    while read sid cram crai; do
      ((k++))
      ((j++))
      # Only print to stdout every 50 samples to limit terminal verbosity
      # This is one of the possible causes of large data egress warnings
      if [ $j -eq 50 ]; then
        echo -e "Checking $wflow for $sid; sample $k of $n $cancer patients"
        j=0
      fi
      ~/code/scripts/check_aou_gatk_status.py \
        --sample-id "$sid" \
        --mode $wflow \
        --bucket $WORKSPACE_BUCKET \
        --status-tsv $status_tsv \
        --update-status \
        --unsafe \
      2> /dev/null
    done < ~/data/cram_paths/$cancer.cram_paths.tsv
  fi
}


# Update status table of sample progress
update_status_table() {
  if [ $# -ne 1 ]; then
    echo "Must specify (gatk-hc|gatk-sv|gvcf-pp) as first positional argument"
  else
    wflow=$1
  fi
  while read cancer; do
    awk -v cancer=$cancer -v wflow=$wflow \
      '{ if ($1!=cancer || $2!=wflow) print $0 }' \
      ~/cromshell/progress/gatk.sample_progress.summary.tsv \
    > ~/cromshell/progress/gatk.sample_progress.summary.tsv2
    cut -f2 \
      ~/cromshell/progress/$cancer.$( echo $wflow | sed 's/-/_/g' ).sample_progress.tsv \
    | sort | uniq -c | sort -nrk1,1 \
    | awk -v cancer=$cancer -v wflow=$wflow -v OFS="\t" \
      '{ print cancer, wflow, $2, $1 }' \
    >> ~/cromshell/progress/gatk.sample_progress.summary.tsv2
    mv \
      ~/cromshell/progress/gatk.sample_progress.summary.tsv2 \
      ~/cromshell/progress/gatk.sample_progress.summary.tsv
  done < ~/cancers.list
  sort -Vk1,3 ~/cromshell/progress/gatk.sample_progress.summary.tsv
}


# Clean up intermediate files generated by Cromwell and collected by check_status
cleanup_garbage() {
  dt_fmt=$( date '+%m_%d_%Y.%Hh%Mm%Ss' )
  garbage_uri="$WORKSPACE_BUCKET/dumpster/dfci_g2c.aou_rw.$dt_fmt.garbage"
  if [ -e ~/uris_to_delete.list ]; then
    gsutil -m cp ~/uris_to_delete.list $garbage_uri
    rm ~/uris_to_delete.list
    echo -e "{\"DeleteGcpObjects.uri_list\": \"$garbage_uri\"}" \
    > ~/cromshell/inputs/empty_dumpster.$dt_fmt.inputs.json
    cromshell --no_turtle -t 120 -mc submit \
      --options-json ~/code/refs/json/aou.cromwell_options.default.json \
      ~/code/wdl/pancan_germline_wgs/DeleteGcpObjects.wdl \
      ~/cromshell/inputs/empty_dumpster.$dt_fmt.inputs.json \
    | jq .id | tr -d '"' \
    >> ~/cromshell/job_ids/empty_dumpster.job_ids.list
  fi
}

