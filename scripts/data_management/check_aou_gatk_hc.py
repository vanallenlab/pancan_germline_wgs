#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Check completion status of GATK-HC for a single All of Us sample
"""


import argparse
import pandas as pd
import subprocess
from os import getenv, path
from re import sub
from sys import stdout


# Global formats
gvcf_fmt = '{}/cromwell/outputs/HaplotypeCallerGvcf_GATK4/{}/call-MergeGVCFs/wgs_{}.g.{}'



def find_gvcfs(workflow_ids, bucket, sid, suffix='vcf.gz'):
    """
    Attempts to find output gVCFs for a list of workflow_ids
    """

    # If no workflow IDs are provided, return empty list
    # This is simply to avoid nested if statements later in the script
    if len(workflow_ids) == 0:
        return []

    # Find GCP URIs for all final output gVCFs
    expected_gvcfs = {gvcf_fmt.format(bucket, wid, sid, suffix) : wid \
                      for wid in workflow_ids}
    gvcf_query = 'gsutil -m ls ' + ' '.join(expected_gvcfs.keys())
    gvcf_query_res = subprocess.run(gvcf_query, capture_output=True, 
                                    shell=True, check=False, text=True)
    gvcfs_found = gvcf_query_res.stdout.rstrip().split('\n')

    # Extract workflow IDs for gVCFs found, if any
    return [uri.split('/')[6] for uri in gvcfs_found if uri.startswith('gs://')]


def check_workflow_status(workflow_id):
    """
    Ping Cromwell server to check status of a single workflow
    """

    cromshell_query = 'cromshell-alpha status ' + workflow_id
    cromshell_query_res = subprocess.run(cromshell_query, capture_output=True, 
                                         shell=True, check=False, text=True)
    return [sub('[,"]', '', s.split('":"')[1]) 
            for s in cromshell_query_res.stdout.split('\n') 
            if 'status' in s][0].lower()


def relocate_outputs(workflow_id, bucket, sid, action='cp'):
    """
    Relocate final gVCF and index for a sample to permanent storage location
    """

    # Relocate gVCF
    src_gvcf = gvcf_fmt.format(bucket, workflow_id, sid, 'vcf.gz')
    dest_gvcf = '{}/dfci-g2c-inputs/aou/gatk-hc/{}.g.vcf.gz'
    subprocess.run(['gsutil', '-m', action, src_gvcf, dest_gvcf], shell=True)

    # Relocate tabix index
    src_tbi = src_gvcf + '.tbi'
    dest_tbi = dest_gvcf + '.tbi'
    subprocess.run(['gsutil', '-m', action, src_tbi, dest_tbi], shell=True)


def collect_trash(workflow_ids, bucket, dumpster_path):
    """
    Add all execution files associated with execution of workflow_ids to the list
    of files to be deleted
    """

    # Find list of all files present in execution buckets
    ex_fmt = '{}/cromwell/execution/HaplotypeCallerGvcf_GATK4/{}/**'
    ex_uris = [ex_fmt.format(bucket, wid) for wid in workflow_ids]
    garbage_raw = subprocess.run('gsutil -m ls ' + ' '.join(ex_uris), 
                                 check=False, capture_output=True, 
                                 shell=True, text=True).stdout
    with open(dumpster_path, 'a') as fout:
        for uri in garbage_raw.rstrip().split('\n'):
            if uri.startswith('gs://'):
                fout.write(uri + '\n')


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-s', '--sample-id', help='Sample ID', required=True)
    parser.add_argument('-b', '--bucket', help='Root bucket [defaut: use ' +
                        '$WORKSPACE_BUCKET environment variable]')
    parser.add_argument('-t', '--status-tsv', help='Two-column .tsv of sample ' +
                        'ID and last known sample status')
    parser.add_argument('-u', '--update-status', default=False, action='store_true',
                        help='Update --status-tsv instead of printing status to stdout')
    parser.add_argument('-r', '--base-directory', default='/home/jupyter',
                        help='Base directory for exectuion. Used to satisfy ' +
                        'assumptions about the locations of Cromshell log files ' + 
                        'and other metadata.')
    parser.add_argument('-d', '--dumpster', default='uris_to_delete.list',
                        help='Path to file collecting all GCP URIs to be deleted. ' +
                        'Note that this script will only ever append to this file, ' +
                        'not overwrite it.')
    args = parser.parse_args()

    # Check to make sure --bucket is set
    bucket = args.bucket
    if args.bucket is None:
        bucket = getenv('WORKSPACE_BUCKET')
    if bucket is None:
        exit('Must provide a valid value for --bucket. See --help for more information.')

    # Default status is always "not_started"
    sid = args.sample_id
    status = 'not_started'

    # If --status-tsv is provided, read last known sample status and use that 
    # information to  speed up the downstream checks
    if args.status_tsv is not None:
        status_tsv = pd.read_csv(args.status_tsv, sep='\t', header=None, 
                                 names='sample_id status'.split())
        status_tsv = status_tsv.astype(str)
        tsv_hit = status_tsv.sample_id == sid
        if tsv_hit.sum() == 1:
            status = status_tsv.loc[tsv_hit, 'status'].values[0]
    else:
        status_tsv = None

    # Check possible statuses in a specific order to minimize necessary 
    # gsutil/cromshell queries
    while True:

        # # First, check if sample output has been already staged
        # # If so, nothing more needs to be done
        # if status == 'staged':
        #     break

        # Second, check if sample has any previous workflow submissions
        # If not, sample has not yet been started and can be reported as such
        workflow_ids_fmt = '{}/cromshell/job_ids/{}.gatk_hc.job_ids.list'
        workflow_ids_path = workflow_ids_fmt.format(args.base_directory, sid)
        if status == 'not_started' \
        and not path.isfile(workflow_ids_path):
            break

        # If sample isn't either unstarted or fully staged, proceed with working 
        # through the possible intermediate states

        # Read ordered list of workflow IDs
        # Order matters because we want to keep the *most recent* complete workflow
        with open(workflow_ids_path) as fin:
            wids = [line.rstrip() for line in fin.readlines()]
            wid_priority = {wid : i + 1 for i, wid in enumerate(wids[::-1])}
        
        # Try to find workflows with final output gVCFs and indexes
        wids_with_gvcf = find_gvcfs(wids, bucket, sid)
        wids_with_tbi = find_gvcfs(wids_with_gvcf, bucket, sid, 'vcf.gz.tbi')

        # If at least one copy of gVCF + index are found, check cromwell status
        # of most recent job ID to confirm job was successful
        for wid in sorted(wids_with_tbi, key=lambda wid: wid_priority[wid]):

            workflow_status = check_workflow_status(wid)

            # If job was successful but gVCF had not been staged yet, relocate 
            # gVCF and index to output bucket and mark all temporary files 
            # for deletion
            if workflow_status == 'succeeded':
                relocate_outputs(wid, bucket, sid)
                collect_trash(wids, bucket, args.dumpster)
                status = 'staged'
                break

        # If at least one workflow was successful, no need to check others
        if status == 'staged':
            break

        # If sample has not been staged, assume it is either running or failed
        # This can be determined by simply checking the status of the 
        # highest priority workflow
        status = check_workflow_status(wids[-1])
        break

    # Report status depending on value of --update-status
    if args.update_status and status_tsv is not None:
        if (status_tsv.sample_id == sid).sum() > 0:
            status_tsv.loc[status_tsv.sample_id == sid, 'status'] = status
        else:
            status_tsv.loc[status_tsv.shape[0], :] = [sid, status]
        status_tsv.to_csv(args.status_tsv, sep='\t', na_rep='NA', 
                          header=False, index=False)
    else:
        stdout.write('{}\t{}\n'.format(sid, status))


if __name__ == '__main__':
    main()

