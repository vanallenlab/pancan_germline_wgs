#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Check completion status of a variant calling workflow for a single All of Us sample
"""


# Import libraries
import argparse
import pandas as pd
import subprocess
from os import getenv, path
from re import sub
from sys import stdout, stderr


# Set global constants
# Note that these definitions make strong assumptions about the structure of 
# WDL/Cromwell execution and output buckets
wdl_names = {'gatk-hc' : 'HaplotypeCallerGvcf_GATK4',
             'gatk-sv' : 'GatherSampleEvidence',
             'gvcf-pp' : 'PostprocessGvcf',
             'read-metrics' : 'CalcReadPairProperties'}
hc_fmts = {'gvcf' : '{}/cromwell/outputs/' + wdl_names['gatk-hc'] + '/{}/call-MergeGVCFs/**wgs_{}.g.{}',
           'dest' : '{}/dfci-g2c-inputs/aou/gatk-hc/{}.g.vcf.gz'}
cov_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-CollectCounts/**{}.counts.tsv.gz',
            'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/coverage/{}.counts.tsv.gz'}
cov_metrics_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-GatherSampleEvidenceMetrics/GatherSampleEvidenceMetrics/*/call-CountsMetrics/**{}.raw-counts.tsv',
                    'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/metrics/{}.raw-counts.tsv'}
pe_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-CollectSVEvidence/CollectSVEvidence/*/call-RunCollectSVEvidence/**{}.pe.txt.gz',
           'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/pesr/{}.pe.txt.gz'}
pe_metrics_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-GatherSampleEvidenceMetrics/GatherSampleEvidenceMetrics/*/call-PEMetrics/**{}.pe-file.tsv',
                   'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/metrics/{}.pe-file.tsv'}
sr_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-CollectSVEvidence/CollectSVEvidence/*/call-RunCollectSVEvidence/**{}.sr.txt.gz',
           'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/pesr/{}.sr.txt.gz'}
sr_metrics_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-GatherSampleEvidenceMetrics/GatherSampleEvidenceMetrics/*/call-SRMetrics/**{}.sr-file.tsv',
                   'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/metrics/{}.sr-file.tsv'}
sd_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-CollectSVEvidence/CollectSVEvidence/*/call-RunCollectSVEvidence/**{}.sd.txt.gz',
           'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/pesr/{}.sd.txt.gz'}
manta_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-Manta/Manta/*/call-RunManta/**{}.manta.vcf.gz',
              'dest' : '{}/dfci-g2c-inputs/aou/manta/{}.manta.vcf.gz'}
manta_metrics_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-GatherSampleEvidenceMetrics/GatherSampleEvidenceMetrics/*/call-Manta_Metrics/manta_**{}.vcf.tsv',
                      'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/metrics/manta_{}.vcf.tsv'}
melt_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-MELT/MELT/*/call-RunMELT/**{}.melt.vcf.gz',
             'dest' : '{}/dfci-g2c-inputs/aou/melt/{}.melt.vcf.gz'}
melt_metrics_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-GatherSampleEvidenceMetrics/GatherSampleEvidenceMetrics/*/call-Melt_Metrics/melt_**{}.vcf.tsv',
                     'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/metrics/melt_{}.vcf.tsv'}
wham_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-Whamg/Whamg/*/call-RunWhamgIncludelist/**{}.wham.vcf.gz',
             'dest' : '{}/dfci-g2c-inputs/aou/wham/{}.wham.vcf.gz'}
wham_metrics_fmts = {'src' : '{}/cromwell/outputs/' + wdl_names['gatk-sv'] + '/{}/call-GatherSampleEvidenceMetrics/GatherSampleEvidenceMetrics/*/call-Wham_Metrics/wham_**{}.vcf.tsv',
                     'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/metrics/wham_{}.vcf.tsv'}
sv_fmts = {'cov' : cov_fmts,
           'cov-metrics' : cov_metrics_fmts,
           'pe' : pe_fmts,
           'pe-metrics' : pe_metrics_fmts,
           'sr' : sr_fmts,
           'sr-metrics' : sr_metrics_fmts,
           'sd' : sd_fmts,
           'manta' : manta_fmts,
           'manta-metrics' : manta_metrics_fmts,
           'melt' : melt_fmts,
           'melt-metrics' : melt_metrics_fmts,
           'wham' : wham_fmts,
           'wham-metrics' : wham_metrics_fmts}
pp_fmts = {'gvcf' : '{}/cromwell/outputs/' + wdl_names['gvcf-pp'] + '/{}/call-Step2/**{}.reblocked.g.{}',
           'dest' : '{}/dfci-g2c-inputs/aou/gatk-hc/reblocked/{}.reblocked.g.vcf.gz'}
read_fmts = {'dest' : '{}/dfci-g2c-inputs/aou/gatk-sv/metrics/{}.read_metrics.tsv'}
formats = {'gatk-hc' : hc_fmts, 
           'gatk-sv' : sv_fmts, 
           'gvcf-pp' : pp_fmts,
           'read-metrics' : read_fmts}
sv_has_index = {'cov' : False,
                'cov-metrics' : False,
                'pe' : True,
                'pe-metrics' : False,
                'sr' : True,
                'sr-metrics' : False,
                'sd' : True,
                'manta' : True,
                'manta-metrics' : False,
                'melt' : True,
                'melt-metrics' : False,
                'wham' : True,
                'wham-metrics' : False}
sv_required = {'cov' : True,
               'cov-metrics' : False,
               'pe' : True,
               'pe-metrics' : False,
               'sr' : True,
               'sr-metrics' : False,
               'sd' : True,
               'manta' : True,
               'manta-metrics' : False,
               'melt' : True,
               'melt-metrics' : False,
               'wham' : True,
               'wham-metrics' : False}


def check_if_staged(bucket, sid, mode, metrics_optional=False):
    """
    Returns true if a sample's outputs are found in the staging bucket
    """

    # List of expected outputs depends on mode
    if mode == 'gatk-sv':
        uris = []
        if metrics_optional:
            required_files = {k : v for k, v in sv_has_index.items() \
                              if sv_required.get(k, False)}
        else:
            required_files = sv_has_index
        for key, has_index in required_files.items():
            uri = formats[mode][key]['dest'].format(bucket, sid)
            uris.append(uri)
            if has_index:
                uris.append(uri + '.tbi')

    else:
        uris = [formats[mode]['dest'].format(bucket, sid)]
        if mode in 'gatk-hc gvcf-pp'.split():
            tbi_uri = uri + '.tbi'
            uris.append(tbi_uri)

    # Check for the presence of all expected URIs
    query = 'gsutil -m ls ' + ' '.join(uris)
    query_res = subprocess.run(query, capture_output=True, shell=True, 
                               check=False, text=True)
    uris_found = query_res.stdout.rstrip().split('\n')

    # Return False unless every output is found in expected location
    return len(set(uris).intersection(set(uris_found))) == len(uris)


def find_gvcfs(workflow_ids, bucket, sid, uri_fmt=formats['gatk-hc']['gvcf'],
               suffix='vcf.gz'):
    """
    Attempts to find output gVCFs for a list of workflow_ids
    """

    # If no workflow IDs are provided, return empty list
    # This is simply to avoid nested if statements later in the script
    if len(workflow_ids) == 0:
        return []

    # Find GCP URIs for all final output gVCFs
    expected_gvcfs = {uri_fmt.format(bucket, wid, sid, suffix) : wid \
                      for wid in workflow_ids}
    gvcf_query = 'gsutil -m ls ' + ' '.join(expected_gvcfs.keys())
    gvcf_query_res = subprocess.run(gvcf_query, capture_output=True, 
                                    shell=True, check=False, text=True)
    gvcfs_found = gvcf_query_res.stdout.rstrip().split('\n')

    # Extract workflow IDs for gVCFs found, if any
    return [uri.split('/')[6] for uri in gvcfs_found if uri.startswith('gs://')]


def find_sv_files(workflow_ids, bucket, sid, uri_fmt):
    """
    Attempts to find one specific SV output file for a list of workflow_ids
    """

    # If no workflow IDs are provided, return empty list
    # This is simply to avoid nested if statements later in the script
    if len(workflow_ids) == 0:
        return []

    # Find GCP URIs for all final output files
    expected_uris = {uri_fmt.format(bucket, wid, sid) : wid \
                     for wid in workflow_ids}
    query = 'gsutil -m ls ' + ' '.join(expected_uris.keys())
    query_res = subprocess.run(query, capture_output=True, shell=True, 
                               check=False, text=True)
    uris_found = query_res.stdout.rstrip().split('\n')

    # Extract workflow IDs for output files found, if any
    return [uri.split('/')[6] for uri in uris_found if uri.startswith('gs://')]


def find_complete_wids(workflow_ids, bucket, sid, mode, metrics_optional=False):
    """
    Find the subset of workflow_ids that have fully complete expected output files
    """

    if mode in 'gatk-hc gvcf-pp'.split():
        wids_with_gvcf = find_gvcfs(workflow_ids, bucket, sid, 
                                    uri_fmt=formats[mode]['gvcf'])
        wids_with_tbi = find_gvcfs(wids_with_gvcf, bucket, sid, 
                                   uri_fmt=formats[mode]['gvcf'],
                                   suffix='vcf.gz.tbi')
        return wids_with_tbi

    elif mode == 'gatk-sv':
        surviving_wids = workflow_ids
        if metrics_optional:
            required_files = {k : v for k, v in sv_has_index.items() \
                              if sv_required.get(k, False)}
        else:
            required_files = sv_has_index
        for key in required_files.keys():
            surviving_wids = find_sv_files(surviving_wids, bucket, sid, 
                                           formats[mode][key]['src'])
        return surviving_wids

    elif mode == 'read-metrics':
        return []
    

def check_workflow_status(workflow_id, max_retries=20, timeout=5):
    """
    Ping Cromwell server to check status of a single workflow
    """

    crom_query_res = ''
    attempts = 0
    while attempts < max_retries:
        crom_query = 'cromshell --no_turtle -t ' + str(timeout) + ' -mc status ' + workflow_id
        crom_query_res = subprocess.run(crom_query, capture_output=True, 
                                        shell=True, check=False, text=True).stdout
        res_splits = [sub('[,"]', '', s.split('":"')[1]) 
                      for s in crom_query_res.split('\n') 
                      if 'status' in s]
        if len(res_splits) > 0:
            return res_splits[0].lower()
        else:
            attempts += 1
    
    if attempts == max_retries:
        msg = 'Failed to get workflow status for {} after {} retries\n'
        stderr.write(msg.format(workflow_id, attempts))
        return 'unknown'        


def relocate_outputs(workflow_id, bucket, staging_bucket, sid, mode, 
                     action='cp', verbose=False):
    """
    Relocate final output files for a sample to permanent storage location
    """

    msg = 'Relocating {} to {}\n'

    if mode in 'gatk-hc gvcf-pp'.split():
        # Relocate gVCF
        src_gvcf = formats[mode]['gvcf'].format(bucket, workflow_id, sid, 'vcf.gz')
        dest_gvcf = formats[mode]['dest'].format(staging_bucket, sid)
        if verbose:
            stdout.write(msg.format(src_gvcf, dest_gvcf))
        subprocess.run(' '.join(['gsutil -m', action, src_gvcf, dest_gvcf]), shell=True)

        # Relocate tabix index
        src_tbi = src_gvcf + '.tbi'
        dest_tbi = dest_gvcf + '.tbi'
        if verbose:
            stdout.write(msg.format(src_tbi, dest_tbi))
        subprocess.run(' '.join(['gsutil -m', action, src_tbi, dest_tbi]), shell=True)

    elif mode == 'gatk-sv':
        for key, has_index in sv_has_index.items():
            src_uri = formats[mode][key]['src'].format(bucket, workflow_id, sid)
            dest_uri = formats[mode][key]['dest'].format(staging_bucket, sid)
            if verbose:
                stdout.write(msg.format(src_uri, dest_uri))
            subprocess.run(' '.join(['gsutil -m', action, src_uri, dest_uri]), shell=True)
            if has_index:
                subprocess.run(' '.join(['gsutil -m', action, src_uri + '.tbi', 
                                         dest_uri + '.tbi']), shell=True)


def collect_trash(workflow_ids, bucket, dumpster_path, mode, sample_id, 
                  staging_bucket):
    """
    Add all execution files associated with execution of workflow_ids to the list
    of files to be deleted
    """

    # Write gsutil-compliant search strings for identifying files to delete
    ex_fmt = '{}/cromwell/execution/{}/{}/**'
    ex_uris = [ex_fmt.format(bucket, wdl_names[mode], wid) for wid in workflow_ids]
    out_fmt = '{}/cromwell/outputs/{}/{}/**'
    out_uris = [out_fmt.format(bucket, wdl_names[mode], wid) for wid in workflow_ids]
    all_uris = ex_uris + out_uris

    # For reblocked gVCFs *only*, add staged raw gVCF to dumpster
    # As of Feb 7, 2024, it was determined that only reblocked gVCFs are necessary
    if mode == 'gvcf-pp':
        raw_gvcf_uri = formats['gatk-hc']['dest'].format(staging_bucket, sample_id)
        all_uris += [raw_gvcf_uri, raw_gvcf_uri + '.tbi']

    # Find list of all files present in execution or output buckets
    # and mark those files for deletion by writing their URIs to the dumpster
    garbage_raw = subprocess.run('gsutil -m ls ' + ' '.join(all_uris), 
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
    parser.add_argument('-m', '--mode', required=True, 
                        help='Check status for which GATK job?', 
                        choices='gatk-hc gatk-sv gvcf-pp read-metrics'.split())
    parser.add_argument('-b', '--bucket', help='Root bucket [defaut: use ' +
                        '$WORKSPACE_BUCKET environment variable]')
    parser.add_argument('-o', '--staging-bucket', help='G2C output staging ' + 
                        'bucket [defaut: same value as --bucket]')
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
    parser.add_argument('--unsafe', default=False, action='store_true',
                        help='Don\'t bother checking Cromwell workflow status ' +
                        'if expected outputs are found. Not recommended.')
    parser.add_argument('--always-print-status-to-stdout', default=False, 
                        action='store_true', help='Always print sample status ' +
                        'to stdout even if --update-status is specified')
    parser.add_argument('--metrics-optional', default=False, action='store_true',
                        help='Treat output metrics files as optional for ' +
                        'determining whether a sample is complete.')
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

        # First, check if sample output has been already staged
        # If so, nothing more needs to be done
        if status == 'staged':
            break
        if check_if_staged(staging_bucket, sid, args.mode, args.metrics_optional):
            status = 'staged'
            break                

        # Second, check if sample has any previous workflow submissions
        # If not, sample has not yet been started and can be reported as such
        workflow_ids_fmt = '{}/cromshell/job_ids/{}.{}.job_ids.list'
        workflow_ids_path = workflow_ids_fmt.format(args.base_directory, sid, 
                                                    sub('-', '_', args.mode))
        if status in 'not_started unknown'.split() \
        and not path.isfile(workflow_ids_path):
            status = 'not_started'
            break

        # If sample isn't either unstarted or fully staged, proceed with working 
        # through the possible intermediate states

        # Read ordered list of workflow IDs
        # Order matters because we want to keep the *most recent* complete workflow
        with open(workflow_ids_path) as fin:
            wids = [line.rstrip() for line in fin.readlines()]
            wid_priority = {wid : i + 1 for i, wid in enumerate(wids[::-1])}

        # For all workflows except for read-metrics, we need to check the
        # execution buckets for complete outputs and stage/cleanup
        if args.mode != 'read-metrics':
            
            # Try to find workflows with full set of complete final outputs
            complete_wids = find_complete_wids(wids, bucket, sid, args.mode, 
                                               args.metrics_optional)

            # If at least one copy of gVCF + index are found, check cromwell status
            # of most recent job ID to confirm job was successful
            for wid in sorted(complete_wids, key=lambda wid: wid_priority[wid]):

                # The exception to this is if --unsafe is set, in which case the
                # workflow is assumed to be successful
                if not args.unsafe:
                    workflow_status = check_workflow_status(wid)
                else:
                    workflow_status = 'succeeded'

                # If job was successful but files have not been staged yet, relocate 
                # files to output bucket and mark all temporary files for deletion
                if workflow_status == 'succeeded':
                    relocate_outputs(wid, bucket, staging_bucket, sid, args.mode)
                    collect_trash(wids, bucket, args.dumpster, args.mode, 
                                  args.sample_id, staging_bucket)
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
    if args.always_print_status_to_stdout \
    or not (args.update_status and status_tsv is not None):
        stdout.write('{}\t{}\n'.format(sid, status))


if __name__ == '__main__':
    main()

