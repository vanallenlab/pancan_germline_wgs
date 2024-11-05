#!/usr/bin/env bash

# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

# Helper bash functions for running the GATK-SV cohort mode pipeline 
# on the All of Us Researcher Workbench


# Note: makes use of helper functions in general_bash_utils.sh, which must be
# sourced before this file is used


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
      export SAMPLES=$( collapse_txt staging/$BATCH/$batch.samples.list )
      while read sid oid cohort; do
        if [ "$cohort" == "aou" ]; then
          gs_base="$MAIN_WORKSPACE_BUCKET/dfci-g2c-inputs"
        else
          gs_base="gs://dfci-g2c-inputs"
        fi
        echo "$gs_base/$cohort/gatk-sv/coverage/$oid.counts.tsv.gz"
      done < staging/$BATCH/$batch.sample_info.tsv \
      > $sub_dir/$BATCH.cov.list
      export COUNTS=$( collapse_txt $sub_dir/$BATCH.cov.list )
      ;;

    04)
      wdl="code/wdl/gatk-sv/GatherBatchEvidence.wdl"
      export SAMPLES=$( collapse_txt staging/$BATCH/$batch.samples.list )
      for suf in cov.list; do
        lfile="$sub_dir/${BATCH}$suf"
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
      export COUNTS=$( collapse_txt $sub_dir/$BATCH.cov.list )
      export PESR_DISC=$( collapse_txt $sub_dir/$BATCH.pe.list )
      export PESR_SD=$( collapse_txt $sub_dir/$BATCH.sd.list )
      export PESR_SPLIT=$( collapse_txt $sub_dir/$BATCH.sd.list )
      export MANTA_VCFS=$( collapse_txt $sub_dir/$BATCH.manta.list )
      export MELT_VCFS=$( collapse_txt $sub_dir/$BATCH.melt.list )
      export WHAM_VCFS=$( collapse_txt $sub_dir/$BATCH.wham.list )
      ;;

    *)
      echo "Module number $module_idx not recognized by submit_gatsv_module. Exiting."
      return
      ;;
  
  esac

  # Format input .json
  ~/code/scripts/envsubst.py \
    -i code/refs/json/dfci-g2c.gatk-sv.$sub_name.inputs.template.json \
  | sed 's/\t//g' | paste -s -d\ \
  > cromshell/inputs/$BATCH.$sub_name.inputs.json

  # Submit job and add job ID to list of jobs for this sample
  cmd="cromshell --no_turtle -t 120 -mc submit"
  cmd="$cmd --options-json code/refs/json/aou.cromwell_options.default.json"
  cmd="$cmd --dependencies-zip gatksv.dependencies.zip"
  cmd="$cmd $wdl cromshell/inputs/$BATCH.$sub_name.inputs.json"
  eval $cmd | jq .id | tr -d '"' \
  >> cromshell/job_ids/$BATCH.$sub_name.job_ids.list
}


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

  # Process batch based on reported status
  status=$( awk -v bid="$bid" -v midx="$module_idx" \
              '{ if ($1==bid && $2==midx) print $3 }' \
              $tracker )
  if [ -z $status ]; then
    echo -e "Failed getting status of batch '$bid' for module '$module_idx'. Exiting."
    return
  fi
  case "$status" in
    "not_started|failed|doomed|unknown")
      submit_batch_module "$bid" "$module_idx"
      ;;
  esac
}
