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
      sub_name="03-TrainGCNV"
      gate=1
      max_resub=3
      ;;
    04)
      sub_name="04-GatherBatchEvidence"
      gate=30
      max_resub=15
      ;;
    05)
      sub_name="05-ClusterBatch"
      gate=2
      max_resub=2
      ;;
    05B)
      sub_name="05B-ExcludeClusteredOutliers"
      gate=0
      max_resub=1
      ;;
    05C)
      sub_name="05C-ReclusterBatch"
      gate=2
      max_resub=2
      ;;
    06)
      sub_name="06-GenerateBatchMetrics"
      gate=30
      max_resub=5
      ;;
    07)
      sub_name="07-FilterBatchSites"
      gate=20
      max_resub=2
      ;;
  esac

  # Process batch based on reported status
  status=$( awk -v bid="$bid" -v midx="$module_idx" \
              '{ if ($1==bid && $2==midx) print $3 }' \
              $tracker )
  if [ -z $status ]; then
    echo -e "Failed getting status of batch '$bid' for module '$module_idx'. Exiting."
    return 1
  fi
  jid_list=~/cromshell/job_ids/$bid.$sub_name.job_ids.list
  if [ -e $jid_list ]; then
    n_prev_subs=$( cat $jid_list | wc -l )
  else
    n_prev_subs=0
  fi
  case "$status" in
    not_started|failed|doomed|aborted)
      if [ $n_prev_subs -le $max_resub ]; then
        echo -e "Submitting GATK-SV module $module_idx for batch $bid"
        submit_batch_module "$bid" "$module_idx"
        echo -e "Waiting $gate minutes before continuing to avoid choking Cromwell server..."
        sleep ${gate}m
      else
        echo -e "Batch $bid has hit the maximum of $max_resub submissions for GATK-SV module $module_idx. Not submitting a new workflow for this batch."
      fi
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

  # Customize gate duration for certain short-running workflows
  outer_gate=60
  case $module_idx in
    05|05C)
      outer_gate=30
      ;;
    05B)
      outer_gate=10
      ;;
    07)
      outer_gate=20
  esac

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
      cleanup_garbage
    done < batch_info/dfci-g2c.gatk-sv.batches.w$WN.list
    echo -e "Finished checking status of all batches for GATK-SV module $module_idx"
    echo -e "Status of all batches for all GATK-SV modules:"
    report_gatksv_status $module_idx
    if [ $( _count_remaining $module_idx ) -eq 0 ]; then
      break
    fi
    echo -e "Waiting $outer_gate minutes before checking progress again..."
    sleep ${outer_gate}m
  done

  echo -e "All batches finished for GATK-SV module $module_idx. Ending monitor routine."
}


# Helper function to report progress for a single module for all batches
report_gatksv_status() {
  # Optional first argument: module index
  tracker="cromshell/progress/gatksv.batch_modules.progress.tsv"
  if [ -e $tracker ]; then
    if [ $# -eq 0 ]; then
      cat $tracker
    else
      awk -v query=$1 -v FS="\t" -v OFS="\t" '{ if ($2==query) print }' $tracker
    fi
  fi
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

  # For modules after 03, check to ensure that the previous module has been staged
  case $module_idx in
    05B)
      prev_idx="04"
      ;;
    05C)
      prev_idx="05B"
      ;;
    06)
      prev_idx="05C"
      ;;
    *)
      prev_idx="0$(( $module_idx - 1 ))"    
      ;;
  esac
  prev_module_outputs_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/$prev_idx/$BATCH/$BATCH.gatksv_module_$prev_idx.outputs.json
  if [ $( gsutil ls $prev_module_outputs_json | wc -l ) -lt 1 ]; then
    echo "Module $module_idx requires the outputs from module $prev_idx to be staged, but staged outputs .json was not found for batch $BATCH. Exiting."
    return 126
  fi

  # Some modules require multiple prior outputs; these extra outputs are checked here
  case $module_idx in
    06)
      # Also requires 04 output
      module_04_outputs_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/04/$BATCH/$BATCH.gatksv_module_04.outputs.json
      if [ $( gsutil ls $module_04_outputs_json | wc -l ) -lt 1 ]; then
        echo "Module $module_idx requires the outputs from module 04 to be staged, but staged outputs .json was not found for batch $BATCH. Exiting."
        return 126
      fi
      ;;
    07)
      # Also requires 05C output
      module_05C_outputs_json=$MAIN_WORKSPACE_BUCKET/dfci-g2c-callsets/gatk-sv/module-outputs/05C/$BATCH/$BATCH.gatksv_module_05C.outputs.json
      if [ $( gsutil ls $module_05C_outputs_json | wc -l ) -lt 1 ]; then
        echo "Module $module_idx requires the outputs from module 05C to be staged, but staged outputs .json was not found for batch $BATCH. Exiting."
        return 126
      fi
      ;;
  esac

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
    05)
      module_name="ClusterBatch"
      ;;
    05B)
      module_name="ExcludeClusteredOutliers"
      ;;
    05C)
      module_name="ReclusterBatch"
      ;;
    06)
      module_name="GenerateBatchMetrics"
      ;;
    07)
      module_name="FilterBatchSites"
      ;;
    *)
      echo "Module number $module_idx not recognized by submit_gatsv_module. Exiting."
      return 2
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
  json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.$sub_name.inputs.template.json
  case $module_idx in

    #############
    # MODULE 03 #
    #############
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

    #############
    # MODULE 04 #
    #############
    04)
      wdl="code/wdl/gatk-sv/GatherBatchEvidence.wdl"

      # Build necessary input arrays
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

      # Begin building .json template updates
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
EOF

      # Get smart about resubmissions: try to determine failure mode and 
      # increase memory for just that failing job
      if [ -e cromshell/job_ids/$BATCH.$sub_name.job_ids.list ]; then
        n_prev_subs=$( cat cromshell/job_ids/$BATCH.$sub_name.job_ids.list | wc -l )
      else
        n_prev_subs=0
      fi
      if [ $n_prev_subs -gt 0 ]; then
        last_sub_id=$( tail -n1 cromshell/job_ids/$BATCH.$sub_name.job_ids.list )
        check_cromwell_return_codes \
          $WORKSPACE_BUCKET/cromwell/execution/$module_name/$last_sub_id \
        > $sub_dir/$BATCH.$last_sub_id.fail_rcs.list
  
        # gCNV case mode
        gcnv_ram=$( jq '."GatherBatchEvidence.runtime_attr_case"' \
                      cromshell/inputs/$BATCH.$sub_name.inputs.json \
                    | fgrep mem_gb | awk '{ print $NF }' )
        gcnv_fails=$( fgrep gCNV $sub_dir/$BATCH.$last_sub_id.fail_rcs.list | wc -l )
        if [ $gcnv_fails -gt 0 ] || [ ! -z $gcnv_ram ]; then
          # Get previous cn.MOPS NormalR1 memory
          if [ -z $gcnv_ram ]; then
            gcnv_ram=48
          fi
          if [ $gcnv_fails -gt 0 ]; then
            # Increment prior memory by +4 GB up to max permitted of 64GB
            gcnv_ram=$(( $gcnv_ram + 4 ))
          fi
          if [ $gcnv_ram -gt 32 ]; then
            gcnv_ram=32
            echo -e "WARNING: module $module_idx for $BATCH is requesting more memory for GATK-gCNV but has already been attempted at max allowed value (32GB)"
          fi
          cat << EOF >> $sub_dir/$BATCH.$sub_name.updates.json
    "GatherBatchEvidence.runtime_attr_case" : { "mem_gb" : $gcnv_ram },
EOF
        fi

        # cn.MOPS
        cnmops_ram=$( jq '."GatherBatchEvidence.cnmops_sample3_runtime_attr"' \
                          cromshell/inputs/$BATCH.$sub_name.inputs.json \
                    | fgrep mem_gb | awk '{ print $NF }' )
        cnmops_fails=$( fgrep CNMOPS $sub_dir/$BATCH.$last_sub_id.fail_rcs.list | wc -l )
        if [ $cnmops_fails -gt 0 ] || [ ! -z $cnmops_ram ]; then
          # Get previous cn.MOPS NormalR1 memory
          if [ -z $cnmops_ram ]; then
            cnmops_ram=48
          fi
          if [ $cnmops_fails -gt 0 ]; then
            # Increment prior memory by +4 GB up to max permitted of 64GB
            cnmops_ram=$(( $cnmops_ram + 4 ))
          fi
          if [ $cnmops_ram -gt 64 ]; then
            cnmops_ram=64
            echo -e "WARNING: module $module_idx for $BATCH is requesting more memory for cn.MOPS but has already been attempted at max allowed value (64GB)"
          fi
          cat << EOF >> $sub_dir/$BATCH.$sub_name.updates.json
    "GatherBatchEvidence.cnmops_sample3_runtime_attr" : { "mem_gb" : $cnmops_ram , "boot_disk_gb" : 20, "disk_gb" : 100 },
EOF
        fi
      fi

      # Finish building .json template updates
      cat << EOF >> $sub_dir/$BATCH.$sub_name.updates.json
    "GatherBatchEvidence.samples" : $( collapse_txt staging/$BATCH/$batch.samples.list ),
    "GatherBatchEvidence.wham_vcfs" : $( collapse_txt $sub_dir/$BATCH.wham.list )
}
EOF
      ;;

    #############
    # MODULE 05 #
    #############
    05)
      wdl="code/wdl/gatk-sv/ClusterBatch.wdl"
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "ClusterBatch.del_bed" : $( gsutil cat $prev_module_outputs_json | jq .merged_dels ),
    "ClusterBatch.dup_bed" : $( gsutil cat $prev_module_outputs_json | jq .merged_dups ),
    "ClusterBatch.manta_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_manta_vcf_tar ),
    "ClusterBatch.melt_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_melt_vcf_tar ),
    "ClusterBatch.wham_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_wham_vcf_tar )
}
EOF
      ;;

    ##############
    # MODULE 05B #
    ##############
    05B)
      wdl="code/wdl/pancan_germline_wgs/ExcludeClusteredOutliers.wdl"
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "ExcludeClusteredOutliers.del_bed" : $( gsutil cat $prev_module_outputs_json | jq .merged_dels ),
    "ExcludeClusteredOutliers.dup_bed" : $( gsutil cat $prev_module_outputs_json | jq .merged_dups ),
    "ExcludeClusteredOutliers.manta_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_manta_vcf_tar ),
    "ExcludeClusteredOutliers.melt_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_melt_vcf_tar ),
    "ExcludeClusteredOutliers.wham_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .std_wham_vcf_tar )
}
EOF
      ;;

    ##############
    # MODULE 05C #
    ##############
    05C)
      wdl="code/wdl/gatk-sv/ClusterBatch.wdl"
      json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.05-ClusterBatch.inputs.template.json
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "ClusterBatch.del_bed" : $( gsutil cat $prev_module_outputs_json | jq .del_bed_cleaned ),
    "ClusterBatch.dup_bed" : $( gsutil cat $prev_module_outputs_json | jq .dup_bed_cleaned ),
    "ClusterBatch.manta_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .manta_vcf_tar_cleaned ),
    "ClusterBatch.melt_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .melt_vcf_tar_cleaned ),
    "ClusterBatch.wham_vcf_tar" : $( gsutil cat $prev_module_outputs_json | jq .wham_vcf_tar_cleaned )
}
EOF
      ;;


    #############
    # MODULE 06 #
    #############
    06)
      wdl="code/wdl/gatk-sv/GenerateBatchMetrics.wdl"
      json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.06-GenerateBatchMetrics.inputs.template.json
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "GenerateBatchMetrics.baf_metrics" : $( gsutil cat $module_04_outputs_json | jq .merged_BAF ),
    "GenerateBatchMetrics.coveragefile" : $( gsutil cat $module_04_outputs_json | jq .merged_bincov ),
    "GenerateBatchMetrics.depth_vcf" : $( gsutil cat $prev_module_outputs_json | jq .clustered_depth_vcf ),
    "GenerateBatchMetrics.discfile" : $( gsutil cat $module_04_outputs_json | jq .merged_PE ),
    "GenerateBatchMetrics.manta_vcf" : $( gsutil cat $prev_module_outputs_json | jq .clustered_manta_vcf ),
    "GenerateBatchMetrics.medianfile" : $( gsutil cat $module_04_outputs_json | jq .median_cov ),
    "GenerateBatchMetrics.melt_vcf" : $( gsutil cat $prev_module_outputs_json | jq .clustered_melt_vcf ),
    "GenerateBatchMetrics.splitfile" : $( gsutil cat $module_04_outputs_json | jq .merged_SR ),
    "GenerateBatchMetrics.wham_vcf" : $( gsutil cat $prev_module_outputs_json | jq .clustered_wham_vcf )
}
EOF
      ;;


    #############
    # MODULE 07 #
    #############
    07)
      wdl="code/wdl/gatk-sv/FilterBatchSites.wdl"
      json_input_template=code/refs/json/gatk-sv/dfci-g2c.gatk-sv.07-FilterBatchSites.inputs.template.json
      cat << EOF > $sub_dir/$BATCH.$sub_name.updates.json
{
    "FilterBatchSites.depth_vcf": $( gsutil cat $module_05C_outputs_json | jq .clustered_depth_vcf ),
    "FilterBatchSites.evidence_metrics": $( gsutil cat $prev_module_outputs_json | jq .metrics ),
    "FilterBatchSites.evidence_metrics_common": $( gsutil cat $prev_module_outputs_json | jq .metrics_common ),
    "FilterBatchSites.manta_vcf": $( gsutil cat $module_05C_outputs_json | jq .clustered_manta_vcf ),
    "FilterBatchSites.melt_vcf": $( gsutil cat $module_05C_outputs_json | jq .clustered_melt_vcf ),
    "FilterBatchSites.wham_vcf": $( gsutil cat $module_05C_outputs_json | jq .clustered_wham_vcf )
}
EOF
      ;;

    *)
      echo "Module number $module_idx not recognized by submit_gatsv_module. Exiting."
      return 2
      ;;
  
  esac

  # Format input .json
  ~/code/scripts/envsubst.py \
    -i $json_input_template \
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

