#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Helper bash functions for running the GATK-SV cohort mode pipeline 
# on the All of Us Researcher Workbench


# Note: makes use of helper functions in general_bash_utils.sh, which must be
# sourced before this file is used


# Check the status of a GATK-SV module for a single batch
# Also updates running tracker
# Also manages submission or staging of outputs as needed
check_batch_module() {
   # Check inputs
  if [ $# -ne 2 ]; then
    echo "Must specify batch ID and [03-08,10] as first two positional arguments"
    return
  else
    export bid=$1
    export module_idx=$2
  fi

  # Check to confirm tracker exists; if it doesn't, create it
  tracker="cromshell/progress/gatksv.batch_modules.progress.tsv"
  if ! [ -e $tracker ]; then touch $tracker; fi

  # Update status of most recent submission for batch
  ~/code/scripts/check_gatksv_batch_module_status.py \
    --batch-id "$bid" \
    --module-index "$module_idx" \
    --staging-bucket "$MAIN_WORKSPACE_BUCKET" \
    --status-tsv "$tracker" \
    --update-status

  # Get module-specific parameters
  case $module_idx in
    03)
      gate=1
      ;;
    04)
      gate=40
      ;;
  esac

  # Process batch based on reported status
  status=$( awk -v bid="$bid" -v midx="$module_idx" \
              '{ if ($1==bid && $2==midx) print $3 }' \
              $tracker )
  if [ -z $status ]; then
    echo -e "Failed getting status of batch '$bid' for module '$module_idx'. Exiting."
    return
  fi
  case "$status" in
    not_started|failed|doomed|unknown)
      echo -e "Submitting GATK-SV module $module_idx for batch $bid"
      submit_batch_module "$bid" "$module_idx"
      echo -e "Waiting $gate minutes before continuing to avoid choking Cromwell server..."
      sleep ${gate}m
      ;;
  esac
}


# Loop to manage submissions & progress tracking for a GATK-SV module for all batches
module_submission_routine_all_batches() {
   # Check input
  if [ $# -ne 1 ]; then
    echo "Must specify [03-08,10] as only positional argument"
    return
  else
    export module_idx=$1
  fi

  WN=$( get_workspace_number )
  batches_list="batch_info/dfci-g2c.gatk-sv.batches.w$WN.list"
  tracker=cromshell/progress/gatksv.batch_modules.progress.tsv

  _count_remaining() {
    if [ -e $tracker ]; then
      incomplete=$( awk -v midx=$1 '{ if ($2==midx && $3!="staged") print }' \
                      $tracker | wc -l )
      missing=$( awk -v midx=$1 '{ if ($2==midx) print $1 }' $tracker \
                 | fgrep -wvf - $batches_list | wc -l )
    else
      incomplete=0
      missing=$( cat $batches_list | wc -l )
    fi
    echo $(( $incomplete + $missing ))
  }

  while [ $( _count_remaining $module_idx ) -gt 0 ]; do
    while read bid; do
      echo -e "Checking status of $bid for GATK-SV module $module_idx"
      check_batch_module $bid $module_idx
    done < batch_info/dfci-g2c.gatk-sv.batches.w$WN.list
    echo -e "Finished checking status of all batches for GATK-SV module $module_idx"
    echo -e "Status of all batches for all GATK-SV modules:"
    cat $tracker
    echo -e "Cleaning up unnecessary Cromwell execution & output files"
    cleanup_garbage
    if [ $( _count_remaining $module_idx ) -eq 0 ]; then
      break
    fi
    echo -e "Waiting 60 minutes before checking progress..."
    sleep 60m
  done

  echo -e "All batches finished for GATK-SV module $module_idx. Ending monitor routine."
}


# Submit a single GATK-SV module for a single batch
submit_batch_module() {
   # Check inputs
  if [ $# -ne 2 ]; then
    echo "Must specify batch ID and [03-08,10] as first two positional arguments"
    return
  else
    export BATCH=$1
    export module_idx=$2
  fi

  # Set constants
  g2c_mfst="dfci-g2c.intake_qc.all.post_qc_batching.tsv.gz"
  batch_sid_list="batch_info/sample_lists/$BATCH.samples.list"

  # Make batch staging directory, if necessary
  for dir in staging "staging/$BATCH"; do
    if ! [ -e $dir ]; then
      mkdir ~/$dir
    fi
  done
  case $module_idx in
    03)
      module_name="TrainGCNV"
      ;;
    04)
      module_name="GatherBatchEvidence"
      prev_module_outputs_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/03/$BATCH/$BATCH.gatksv_module_03.outputs.json
      if [ $( gsutil ls $prev_module_outputs_json | wc -l ) -lt 1 ]; then
        echo "Module 04 requires the outputs from module 03 to be staged, but staged outputs .json was not found for batch $BATCH. Exiting."
      fi
      ;;
    *)
      echo "Module number $module_idx not recognized by submit_gatsv_module. Exiting."
      return
      ;;
  esac
  sub_name="${module_idx}-$module_name"
  sub_dir="staging/$BATCH/$sub_name"
  if ! [ -e $sub_dir ]; then mkdir $sub_dir; fi

  # Get map of G2C ID, original ID, and original cohort for this batch
  zcat $g2c_mfst | fgrep -wf $batch_sid_list | cut -f1-3 \
  > staging/$BATCH/$batch.sample_info.tsv
  cut -f1 staging/$BATCH/$batch.sample_info.tsv \
  > staging/$BATCH/$batch.samples.list

  # Set workflow-specific parameters
  case $module_idx in

    03)
      wdl="code/wdl/gatk-sv/TrainGCNV.wdl"
      while read sid oid cohort; do
        if [ "$cohort" == "aou" ]; then
          gs_base="$MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs"
        else
          gs_base="gs://dfci-g2c-inputs"
        fi
        echo "$gs_base/$cohort/gatk-sv/coverage/$oid.counts.tsv.gz"
      done < staging/$BATCH/$batch.sample_info.tsv \
      > $sub_dir/$BATCH.cov.list
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "GatherBatchEvidence.counts" : $( collapse_txt $sub_dir/$BATCH.cov.list ),
    "GatherBatchEvidence.samples" : $( collapse_txt staging/$BATCH/$batch.samples.list )
}
EOF
      ;;

    04)
      wdl="code/wdl/gatk-sv/GatherBatchEvidence.wdl"
      for suf in cov.list pe.list sr.list sd.list manta.list melt.list wham.list; do
        lfile="$sub_dir/${BATCH}.$suf"
        if [ -e $lfile ]; then rm $lfile; fi
      done
      while read sid oid cohort; do
        if [ "$cohort" == "aou" ]; then
          gs_base="$MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs"
        else
          gs_base="gs://dfci-g2c-inputs"
        fi
        echo "$gs_base/$cohort/gatk-sv/coverage/$oid.counts.tsv.gz" >> $sub_dir/$BATCH.cov.list
        echo "$gs_base/$cohort/gatk-sv/pesr/$oid.pe.txt.gz" >> $sub_dir/$BATCH.pe.list
        echo "$gs_base/$cohort/gatk-sv/pesr/$oid.sr.txt.gz" >> $sub_dir/$BATCH.sr.list
        echo "$gs_base/$cohort/gatk-sv/pesr/$oid.sd.txt.gz" >> $sub_dir/$BATCH.sd.list
        echo "$gs_base/$cohort/manta/$oid.manta.vcf.gz" >> $sub_dir/$BATCH.manta.list
        echo "$gs_base/$cohort/melt/$oid.melt.vcf.gz" >> $sub_dir/$BATCH.melt.list
        echo "$gs_base/$cohort/wham/$oid.wham.vcf.gz" >> $sub_dir/$BATCH.wham.list
      done < staging/$BATCH/$batch.sample_info.tsv
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "GatherBatchEvidence.PE_files" : $( collapse_txt $sub_dir/$BATCH.pe.list ),
    "GatherBatchEvidence.SD_files" : $( collapse_txt $sub_dir/$BATCH.sd.list ),
    "GatherBatchEvidence.SR_files" : $( collapse_txt $sub_dir/$BATCH.sr.list ),
    "GatherBatchEvidence.contig_ploidy_model_tar" : $( gsutil cat $prev_module_outputs_json \
                                                       | jq .cohort_contig_ploidy_model_tar \
                                                       | awk '{ if ($1 ~ /gs/) print $1 }' ),
    "GatherBatchEvidence.counts" : $( collapse_txt $sub_dir/$BATCH.cov.list ),
    "GatherBatchEvidence.gcnv_model_tars" : $( gsutil cat $prev_module_outputs_json \
                                               | jq .cohort_gcnv_model_tars \
                                               | paste -s - -d\ | tr -s ' ' ),
    "GatherBatchEvidence.manta_vcfs" : $( collapse_txt $sub_dir/$BATCH.manta.list ),
    "GatherBatchEvidence.melt_vcfs" : $( collapse_txt $sub_dir/$BATCH.melt.list ),
    "GatherBatchEvidence.samples" : $( collapse_txt staging/$BATCH/$batch.samples.list ),
    "GatherBatchEvidence.wham_vcfs" : $( collapse_txt $sub_dir/$BATCH.wham.list )
}
EOF
      ;;

    *)
      echo "Module number $module_idx not recognized by submit_gatsv_module. Exiting."
      return
      ;;
  
  esac

  # Format input .json
  ~/code/scripts/envsubst.py \
    -i code/refs/json/gatk-sv/dfci-g2c.gatk-sv.$sub_name.inputs.template.json \
  | ~/code/scripts/update_json.py \
    -i stdin \
    -u $sub_dir/$BATCH.$sub_name.updates.json \
    -o cromshell/inputs/$BATCH.$sub_name.inputs.json

  # Submit job and add job ID to list of jobs for this sample
  cmd="cromshell --no_turtle -t 120 -mc submit"
  cmd="$cmd --options-json code/refs/json/aou.cromwell_options.default.json"
  cmd="$cmd --dependencies-zip gatksv.dependencies.zip"
  cmd="$cmd $wdl cromshell/inputs/$BATCH.$sub_name.inputs.json"
  eval $cmd | jq .id | tr -d '"' \
  >> cromshell/job_ids/$BATCH.$sub_name.job_ids.list
}

