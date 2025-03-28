#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Manage submission of a single GATK-SV workflow for a single G2C batch
"""


# Import libraries
import argparse
import g2cpy
import json
import pandas as pd
import subprocess
from os import getenv, path, remove
from re import sub
from sys import stdout, stderr


# Set global constants
# Note that these definitions make strong assumptions about the structure of 
# WDL/Cromwell execution and output buckets
wdl_names = {'03' : 'TrainGCNV',
             '04' : 'GatherBatchEvidence',
             '05' : 'ClusterBatch',
             '05B' : 'ExcludeClusteredOutliers',
             '05C' : 'ClusterBatch',
             '06' : 'GenerateBatchMetrics',
             '07' : 'FilterBatchSites',
             '08' : 'FilterBatchSamples',
             '10' : 'GenotypeBatch',
             '14A' : 'FilterCoverageSamples'}
module_alt_name = {'05C' : '05C-ReclusterBatch'}
output_bucket_fmt = '{0}/dfci-g2c-callsets/gatk-sv/module-outputs/{1}/{2}'
output_json_fname_fmt = '{2}.gatksv_module_{1}.outputs.json'
output_json_fmt = '/'.join([output_bucket_fmt, output_json_fname_fmt])
keep_03_outs = 'cohort_contig_ploidy_model_tar cohort_gcnv_model_tars'.split()
keep_04_outs = 'manta_tloc median_cov merged_BAF merged_BAF_index merged_PE ' + \
               'merged_PE_index merged_SR merged_SR_index merged_bincov ' + \
               'merged_bincov_index merged_dels merged_dups std_manta_vcf_tar ' + \
               'std_melt_vcf_tar std_wham_vcf_tar'
keep_04_outs = keep_04_outs.split()
keep_05_outs = 'clustered_depth_vcf clustered_depth_vcf_index clustered_manta_vcf ' + \
               'clustered_manta_vcf_index clustered_melt_vcf ' + \
               'clustered_melt_vcf_index clustered_wham_vcf ' + \
               'clustered_wham_vcf_index clustered_outlier_samples_with_reason ' + \
               'clustered_sv_counts clustered_sv_count_plots metrics_file_clusterbatch'
keep_05_outs = keep_05_outs.split()
keep_05B_outs = 'del_bed_cleaned dup_bed_cleaned manta_vcf_tar_cleaned ' + \
                'melt_vcf_tar_cleaned wham_vcf_tar_cleaned'
keep_05B_outs = keep_05B_outs.split()
keep_05C_outs = keep_05_outs
keep_06_outs = 'metrics metrics_common metrics_file_batchmetrics'.split()
keep_07_outs = 'sites_filtered_melt_vcf sites_filtered_wham_vcf cutoffs ' + \
               'sites_filtered_sv_counts sites_filtered_manta_vcf ' + \
               'sites_filtered_depth_vcf scores'
keep_07_outs = keep_07_outs.split()
keep_08_outs = 'outlier_filtered_depth_vcf outlier_filtered_depth_vcf_index ' + \
               'outlier_filtered_pesr_vcf outlier_filtered_pesr_vcf_index ' + \
               'filtered_batch_samples_file'
keep_08_outs = keep_08_outs.split()
keep_10_outs = 'trained_genotype_pesr_depth_sepcutoff sr_background_fail ' + \
               'trained_SR_metrics genotyped_depth_vcf_index ' + \
               'trained_genotype_depth_pesr_sepcutoff genotyped_pesr_vcf_index ' + \
               'sr_bothside_pass genotyped_depth_vcf metrics_file_genotypebatch ' + \
               'trained_genotype_depth_depth_sepcutoff genotyped_pesr_vcf ' + \
               'regeno_coverage_medians trained_genotype_pesr_pesr_sepcutoff ' + \
               'trained_PE_metrics'
keep_10_outs = keep_10_outs.split()
keep_14A_outs = 'filtered_bincov filtered_bincov_idx filtered_median_cov'.split()
keep_output_keys = {'03' : keep_03_outs,
                    '04' : keep_04_outs,
                    '05' : keep_05_outs,
                    '05B' : keep_05B_outs,
                    '05C' : keep_05C_outs,
                    '06' : keep_06_outs,
                    '07' : keep_07_outs,
                    '08' : keep_08_outs,
                    '10' : keep_10_outs,
                    '14A' : keep_14A_outs}


def check_if_staged(bucket, bid, module_index):
    """
    Returns true if an outputs .json file is found in the expected location

    Note that we bookkeep GATK-SV module outputs using a reserved .json file,
    which should *only* be present if the module was successful and all of the 
    module outputs were able to be staged in --staging-bucket.
    """

    return g2cpy.check_gcp_uris([output_json_fmt.format(bucket, module_index, bid)])


def relocate_outputs(workflow_id, staging_bucket, bid, module_index, 
                     action='cp', timeout=30, max_retries=20, verbose=False):
    """
    Relocate final output files for a module to a permanent storage location
    """

    msg = 'Relocating {} to {}\n'

    wdl_name = wdl_names[module_index]

    # Use cromshell list-outputs to get and read a melted dict-like .txt of outputs
    crom_query = 'cromshell --no_turtle -t ' + str(timeout) + \
                 ' --machine_processable list-outputs ' + workflow_id
    attempts = 0
    while attempts < max_retries:
        crom_query_res = subprocess.run(crom_query, capture_output=True, shell=True, 
                                        check=False, text=True).stdout
        if crom_query_res != '':
            break
        else:
            attempts += 1
        if attempts == max_retries:
            msg = 'cromshell list-outputs failed after {:,} retries with {:,}s timeout'
            exit(msg.format(max_retries, timeout))

    # Iterate over each output mentioned by cromwell and move each to staging_bucket
    # While relocating output files, build an outputs .json to reflect their new locations
    outputs_dict = {}
    keep_keys = keep_output_keys[module_index]
    for entry in crom_query_res.rstrip().split('\n'):

        key, src_uri = sub('^' + wdl_name + '\.', '', entry).split(': ')
    
        if key not in keep_keys:
            continue
    
        # Format destination URI to keep nested structure of subworkflows
        # but to remove all arbitrary Cromwell hashes
        dest_parts = [sub('^call-', '', x) for x in src_uri.split('/') 
                      if x.startswith('call-') or x.startswith('shard-')]
        dest_prefix = output_bucket_fmt.format(staging_bucket, module_index, bid)
        dest_uri = '/'.join([dest_prefix] + dest_parts + [path.basename(src_uri)])
        
        # Stage output
        g2cpy.relocate_uri(src_uri, dest_uri, verbose=verbose)

        # Update outputs .json
        if key not in outputs_dict.keys():
            outputs_dict[key] = list()
        outputs_dict[key].append(dest_uri)

    # Write updated output .json to file and copy that file into staging_bucket
    for key in outputs_dict.keys():
        if len(outputs_dict[key]) == 1:
            outputs_dict[key] = outputs_dict[key][0]
    json_fout = output_json_fname_fmt.format('', module_index, bid)
    with open(json_fout, 'w') as fout:
        json.dump(outputs_dict, fout)
    json_out_uri = output_json_fmt.format(staging_bucket, module_index, bid)
    g2cpy.relocate_uri(json_fout, json_out_uri)
    remove(json_fout)


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--batch-id', help='Batch ID', required=True)
    parser.add_argument('-m', '--module-index', required=True, help='Module number')
    parser.add_argument('--bucket', help='Root bucket [defaut: use ' +
                        '$WORKSPACE_BUCKET environment variable]')
    parser.add_argument('-o', '--staging-bucket', help='G2C output staging ' + 
                        'bucket [defaut: same value as --bucket]')
    parser.add_argument('-t', '--status-tsv', help='Three-column .tsv of batch ' +
                        'ID, workflow name, and last known status')
    parser.add_argument('-u', '--update-status-tsv', default=False, action='store_true',
                        help='Update --status-tsv instead of printing status to stdout')
    parser.add_argument('-r', '--base-directory', default='/home/jupyter',
                        help='Base directory for exectuion. Used to satisfy ' +
                        'assumptions about the locations of Cromshell log files ' + 
                        'and other metadata.')
    parser.add_argument('-d', '--dumpster', default='uris_to_delete.list',
                        help='Path to file collecting all GCP URIs to be deleted. ' +
                        'Note that this script will only ever append to this file, ' +
                        'not overwrite it.')
    parser.add_argument('--always-print-status-to-stdout', default=False, 
                        action='store_true', help='Always print sample status ' +
                        'to stdout even if --update-status is specified')
    args = parser.parse_args()

    # Check to make sure --bucket is set
    bucket = args.bucket
    if args.bucket is None:
        bucket = getenv('WORKSPACE_BUCKET')
    if bucket is None:
        exit('Must provide a valid value for --bucket. See --help for more information.')
    if args.staging_bucket is not None:
        staging_bucket = args.staging_bucket
    else:
        staging_bucket = bucket

    # Default status is always "not_started"
    bid = args.batch_id
    status = 'not_started'

    # Confirm that module is recognized
    wdl_name = wdl_names.get(args.module_index)
    if wdl_name is None:
        msg = 'Module index {} is not currently supported. Exiting.'
        exit(msg.format(args.module_index))
    module_name = module_alt_name.get(args.module_index, None)
    if module_name is None:
        module_name = '-'.join([args.module_index, wdl_name])

    # If --status-tsv is provided, read last known batch status and use that 
    # information to  speed up the downstream checks
    if args.status_tsv is not None:
        status_tsv = pd.read_csv(args.status_tsv, sep='\t', header=None, dtype=str, 
                                 names='batch_id module status'.split())
        status_tsv = status_tsv.astype(str)
        tsv_hit = (status_tsv.batch_id == bid) & (status_tsv.module == args.module_index)
        if tsv_hit.sum() == 1:
            status = status_tsv.loc[tsv_hit, 'status'].values[0]
    else:
        status_tsv = None

    # Check possible statuses in a specific order to minimize necessary 
    # gsutil/cromshell queries
    while True:

        # First, check if batch output has been already staged
        # If so, nothing more needs to be done
        if status == 'staged':
            break
        if check_if_staged(staging_bucket, bid, args.module_index):
            status = 'staged'
            break

        # Second, check if batch has any previous workflow submissions
        # If not, batch has not yet been started and can be reported as such
        workflow_ids_fmt = '{}/cromshell/job_ids/{}.{}.job_ids.list'
        workflow_ids_path = workflow_ids_fmt.format(args.base_directory, bid, 
                                                    module_name)
        if status in 'not_started unknown'.split() \
        and not path.isfile(workflow_ids_path):
            status = 'not_started'
            break

        # If batch isn't either unstarted or fully staged, proceed with working 
        # through the possible intermediate states

        # Read ordered list of workflow IDs
        with open(workflow_ids_path) as fin:
            wids = [line.rstrip() for line in fin.readlines()]

        # Check to make sure at least one workflow ID is reported
        if len(wids) == 0:
            status = 'not_started'
            break

        # Update status according to most recent workflow
        wid = wids[-1]
        status = g2cpy.check_workflow_status(wid, timeout=120)

        # If most recent workflow was successful, stage outputs and clear all 
        # files from Cromwell execution & output buckets
        if status == 'succeeded':
            relocate_outputs(wid, staging_bucket, bid, args.module_index)
            status = 'staged'
            g2cpy.collect_workflow_trash(wids, bucket, wdl_name, args.dumpster)

        # Exit loop after processing most recent workflow ID
        break

    # Report status depending on value of --update-status
    if args.update_status_tsv and status_tsv is not None:
        if tsv_hit.sum() > 0:
            status_tsv.loc[tsv_hit, 'status'] = status
        else:
            status_tsv.loc[status_tsv.shape[0], :] = [bid, args.module_index, status]
        status_tsv.to_csv(args.status_tsv, sep='\t', na_rep='NA', 
                          header=False, index=False)
    if args.always_print_status_to_stdout \
    or not (args.update_status_tsv and status_tsv is not None):
        stdout.write('\t'.join([bid, args.module_index, status]) + '\n')


if __name__ == '__main__':
    main()

