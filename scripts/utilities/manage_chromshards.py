#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Routine to submit, track, manage, and stage a method parallelized per chromosome
"""


# Import libraries
import argparse
import g2cpy
import json
import pandas as pd
import subprocess
import textwrap
from datetime import datetime
from os import environ, getenv, path, putenv, remove
from re import sub
from sys import stdout, stderr
from time import sleep


# Define constants
hg38_primary_contigs = ['chr{}'.format(x+1) for x in range(22)] + ['chrX', 'chrY']


def clean_date():
    """
    Wrapper for clean date & time reporting
    """

    return datetime.now().strftime("%Y-%m-%d at %H:%M:%S")


def startup_report(inputs):
    """
    Prints a startup message to stdout with all runtime info
    """

    start_log="""

    #=======================================================#
    |      DFCI G2C chromosome-sharded workflow manager     |
    | Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu> |
    #=======================================================#

    Launched on {}

    Inputs:"""
    start_log = start_log.format(clean_date())
    print(textwrap.dedent(start_log))

    for key, val in inputs.items():
        print('  - {} : {}'.format(key, val))
    print('')


def report_status(all_status, title=None):
    """
    Summarizes workflow status across all chromosomes held in all_status tracker
    """

    # Group contigs by status
    status_groups = {}
    for k, status in all_status.items():
        if status not in status_groups:
            status_groups[status] = set()
        status_groups[status].add(k)

    # Report status in decreasing frequency
    if title is None:
        print('[{}] Status summary for all contigs:'.format(clean_date()))
    else:
        print(title)
    for status, contigs in sorted(status_groups.items(), 
                                  key=lambda x: len(x[1]), reverse=True):
        print('  - {:,} {}: {}'.format(len(contigs), status, 
                                       ', '.join(g2cpy.chromsort(list(contigs)))))
    print('')


def relocate_outputs(workflow_id, staging_bucket, wdl_name, output_json_uri,
                     action='cp', timeout=30, max_retries=20, verbose=False):
    """
    Relocate final output files for a single workflow to a permanent storage location
    """

    msg = 'Relocating {} to {}\n'

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
    for entry in crom_query_res.rstrip().split('\n'):

        key, src_uri = sub('^' + wdl_name + '\.', '', entry).split(': ')
    
        # Format destination URI to keep nested structure of subworkflows
        # but to remove all arbitrary Cromwell hashes
        dest_parts = [sub('^call-', '', x) for x in src_uri.split('/') 
                      if x.startswith('call-') or x.startswith('shard-')]
        dest_uri = '/'.join([staging_bucket] + dest_parts + [path.basename(src_uri)])
        
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
    json_fout = path.basename(output_json_uri)
    with open(json_fout, 'w') as fout:
        json.dump(outputs_dict, fout)
    g2cpy.relocate_uri(json_fout, output_json_uri)
    remove(json_fout)


def submit_workflow(contig, wdl, input_template, input_json, prev_wids, quiet=False):
    """
    Submit a workflow for a single contig to Cromwell and store the workflow ID in a log
    """

    # Substitute environment variables from input template to specific input
    for vbl in 'CHR CHROM CONTIG chr chrom contig'.split():
        environ[vbl] = contig
        with open(input_template) as fin:
            with open(input_json, 'w') as fout:
                for line in fin.readlines():
                    fout.write(path.expandvars(line))

    # Build the workflow submission command
    cmd = 'cromshell --no_turtle -t 120 -mc submit'
    cmd += ' --options-json code/refs/json/aou.cromwell_options.default.json'
    cmd += ' --dependencies-zip gatksv.dependencies.zip ' + wdl + ' ' + input_json

    # Submit the workflow
    sub_res = subprocess.run(cmd, capture_output=True, shell=True, 
                             check=False, text=True)
    try:
        new_wid = json.loads(sub_res.stdout)['id']
        with open(prev_wids, 'a') as fout:
            fout.write(new_wid + '\n')
        status = 'submitted'
        if not quiet:
            msg = '[{}] Submitted {} for contig {} (workflow ID: {})'
            print(msg.format(clean_date(), path.basename(wdl), contig, new_wid))
    except:
        status = 'submission_error'

    return status


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-w', '--wdl', help='WDL file to execute', required=True)
    parser.add_argument('-i', '--input-json-template', help='Template for input ' +
                        '.json to pass to Cromwell. All environment variables ' +
                        'present will be automatically expanded. Any of $CHR, ' +
                        '$CHROM, $CONTIG, $chr, $chrom, or $contig are reserverd ' +
                        'and will be replaced by contig name within this routine.',
                        required=True)
    parser.add_argument('-o', '--staging-bucket', help='G2C output staging ' + 
                        'bucket. Must provide terminal bucket within which ' +
                        'each chromosome will have outputs staged in separate ' +
                        'directory.')
    parser.add_argument('-n', '--name', help='Method name or nickname. If ' +
                        'not provided, the basename of --wdl will be used.')
    parser.add_argument('--bucket', help='Root bucket [defaut: use ' +
                        '$WORKSPACE_BUCKET environment variable]')
    parser.add_argument('--contig-list', help='List of contig names to ' + 
                        'evaluate [default: all 24 primary autosomes]')
    parser.add_argument('-t', '--status-tsv', help='Two-column .tsv of contig ' +
                        'name and last known status')
    parser.add_argument('-r', '--base-directory', default='/home/jupyter',
                        help='Base directory for exectuion. Used to satisfy ' +
                        'assumptions about the locations of Cromshell log files ' + 
                        'and other metadata.')
    parser.add_argument('-d', '--dumpster', default='uris_to_delete.list',
                        help='Path to file collecting all GCP URIs to be deleted. ' +
                        'Note that this script will only ever append to this file, ' +
                        'not overwrite it.')
    parser.add_argument('-m', '--max-attempts', type=int, default=3, 
                        help='Maximum number of submission attempts for any one contig')
    parser.add_argument('-g', '--gate', type=float, default=20, help='Number of ' +
                        'minutes to wait between monitor cycles')
    parser.add_argument('--workflow-id-log-prefix', help='Optional prefix for ' +
                        'logger files for submitted workflow IDs')
    parser.add_argument('--quiet', default=False, action='store_true',help='Do ' +
                        'not print logging, progress, and diagnostics to stdout')
    args = parser.parse_args()

    # Infer method name if not provided
    wdl_name = path.splitext(path.basename(args.wdl))[0]
    if args.name is None:
        method_name = wdl_name
    else:
        method_name = args.name

    # Check to make sure --bucket is set
    bucket = args.bucket
    if args.bucket is None:
        bucket = getenv('WORKSPACE_BUCKET')
    if bucket is None:
        exit('Must provide a valid value for --bucket. See --help for more information.')

    # Infer expected paths and other contig-specific constants
    input_json_fmt = sub('[/]+$', '', args.base_directory) + \
                     '/cromshell/inputs/' + \
                     sub('.template', '', 
                         path.splitext(path.basename(args.input_json_template))[0]) + \
                     '.{}.json'
    output_bucket_fmt = sub('[/]+$', '', args.staging_bucket) + '/{0}'
    output_json_fmt = output_bucket_fmt + '/' + method_name + '.{0}.outputs.json'
    if args.workflow_id_log_prefix is not None:
        prev_wid_fname_base = args.workflow_id_log_prefix + '.' + method_name
    else:
        prev_wid_fname_base = method_name
    prev_wids_fmt = sub('[/]+$', '', args.base_directory) + \
                    '/cromshell/job_ids/' + prev_wid_fname_base + \
                    '.{}.job_ids.list'

    # Read custom list of contigs, if optioned
    if args.contig_list is None:
        contigs = hg38_primary_contigs
    else:
        with open(args.contig_list) as clist:
            contigs = [x.rstrip() for x in clist.readlines()]
    contigs = g2cpy.chromsort(contigs)

    # Report startup conditions
    if not args.quiet:
        startup_report({'Method name' : method_name,
                        'WDL' : args.wdl,
                        'Input .json template' : args.input_json_template,
                        'Workspace bucket' : bucket,
                        'Staging bucket' : args.staging_bucket,
                        'Local root' : args.base_directory,
                        'Dumpster' : args.dumpster,
                        'Contigs' : ', '.join(contigs),
                        'Quiet' : args.quiet})

    # If --status-tsv is provided and exists, load this file as a dict
    if args.status_tsv is not None \
    and path.exists(args.status_tsv):
        all_status = pd.read_csv(args.status_tsv, sep='\t', header=None).\
                        set_index(0)[1].to_dict()
    # Otherwise, create a dictionary for tracking status
    else:
        all_status = {k : 'unknown' for k in contigs}

    # Report beginning status
    if not args.quiet:
        report_status(all_status, title='Status at launch:')

    # Loop infinitely while any status is not 'staged'
    while True:
        if not args.quiet:
            msg = '[{}] Not all contigs yet staged. Entering another cycle of ' + \
                  'submission management routine.'
            print(msg.format(clean_date()))

        # Loop over all chromosomes
        for contig in contigs:

            # Get most recent status for this contig
            status = all_status.get(contig, 'unknown')

            # Format expected input + output paths
            input_json = input_json_fmt.format(contig)
            prev_wids = prev_wids_fmt.format(contig)
            output_bucket = output_bucket_fmt.format(contig)
            output_json_uri = output_json_fmt.format(contig)

            # Check possible statuses in a specific order to minimize necessary 
            # gsutil/cromshell queries
            while True:

                # First, check if contig output has been already staged
                # If so, nothing more needs to be done
                if status == 'staged':
                    break
                if g2cpy.check_gcp_uris([output_json_uri]):
                    status = 'staged'
                    break

                # Second, check if contig has any previous workflow submissions
                # If not, contig has not yet been started and can be reported as such
                if status in 'not_started unknown'.split() \
                and not path.isfile(prev_wids):
                    status = 'not_started'
                    n_prior_subs = 0
                    break

                # If contig isn't either unstarted or fully staged, proceed with working 
                # through the possible intermediate states

                # Read ordered list of workflow IDs
                if path.exists(prev_wids):
                    with open(prev_wids) as fin:
                        wids = [line.rstrip() for line in fin.readlines()]
                        n_prior_subs = len(wids)
                else:
                    wids = []
                    n_prior_subs = 0

                # Check to make sure at least one workflow ID is reported
                if n_prior_subs == 0:
                    status = 'not_started'
                    break

                # Update status according to most recent workflow
                wid = wids[-1]
                status = g2cpy.check_workflow_status(wid, timeout=120)

                # If most recent workflow was successful, stage outputs and clear all 
                # files from Cromwell execution & output buckets
                if status == 'succeeded':
                    relocate_outputs(wid, output_bucket, wdl_name, output_json_uri)
                    status = 'staged'
                    g2cpy.collect_workflow_trash(wids, bucket, wdl_name,
                                                 args.dumpster)
                    subprocess.run('. code/refs/general_bash_utils.sh ' + \
                                   '&& cleanup_garbage', shell=True, check=False, 
                                   text=True, executable='/bin/bash')
                    if not args.quiet:
                        print('')

                # Exit loop after processing most recent workflow ID
                break

            # Submit new workflow if not already running or staged
            if status not in 'staged submitted running'.split():
                if n_prior_subs < args.max_attempts:
                    status = submit_workflow(contig, args.wdl, args.input_json_template, 
                                             input_json, prev_wids, args.quiet)
                else:
                    if not args.quiet:
                        msg = '[{}] Contig {} has reached the maximum number of ' + \
                              'attempts ({:,}) for {}; skipping.\n'
                        print(msg.format(clean_date(), contig, 
                                         args.max_attempts, method_name))
                    status = 'exhausted'

            # Update & report status
            all_status[contig] = status
            if not args.quiet:
                msg = '[{}] Status of contig {}: {}'
                print(msg.format(clean_date(), contig, status))
            if args.status_tsv is not None:
                with open(args.status_tsv, 'w') as fout:
                    for k, v in all_status.items():
                        fout.write('\t'.join([k, v]) + '\n')

        # Check at end of loop if all contigs are staged
        if len(set(all_status.values()).difference({'staged'})) == 0:
            break

        # Wait before checking again
        else:
            if not args.quiet:
                report_status(all_status)
                msg = '[{}] Finished checking status for all {:,} contigs. ' + \
                      'Waiting {:,} minutes before checking again...\n'
                print(msg.format(clean_date(), len(contigs), args.gate))
            sleep(60 * args.gate)

    # Report completion once all contigs are staged
    msg = '[{}] All contig outputs have been staged! Exiting management routine.'
    print(msg.format(clean_date()))


if __name__ == '__main__':
    main()

