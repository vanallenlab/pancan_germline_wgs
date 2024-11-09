#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Helper functions for data processing on the All of Us Researcher Workbench
"""


import subprocess
from re import sub
from sys import stdout, stderr


def check_workflow_status(workflow_id, max_retries=20, timeout=30):
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


def check_gcp_uris(uris, return_uri_list=False):
    """
    Check if a list of Google URIs exists
    If return_uri_list is True, will return a list of all URIs found
    Otherwise, will return a simple True/False if *all* URIs were found
    """

    # Build URI query
    query = 'gsutil -m ls ' + ' '.join(uris)
    query_res = subprocess.run(query, capture_output=True, shell=True, 
                               check=False, text=True)
    uris_found = query_res.stdout.rstrip().split('\n')

    # Return False unless every URI is found
    if return_uri_list:
        return uris_found
    else:
        return len(set(uris).intersection(set(uris_found))) == len(uris)


def collect_gcp_garbage(uris, dumpster_path='uris_to_delete.list'):
    """
    Writes a list of URIs to delete to a dumpster (flat text) file
    Note that the input URIs need not exist, and can include wildcards/double wildcards
    """

    garbage_raw = subprocess.run('gsutil -m ls ' + ' '.join(uris), 
                                 check=False, capture_output=True, 
                                 shell=True, text=True).stdout
    with open(dumpster_path, 'a') as fout:
        for uri in garbage_raw.rstrip().split('\n'):
            if uri.startswith('gs://'):
                fout.write(uri + '\n')


def relocate_uri(src_uri, dest_uri, action='cp', verbose=False):
    """
    Moves or copies a GCP object from a source URI to a destination URI
    """

    if action not in 'cp mv'.split():
        msg = 'Action "{}" is not recognized. Please specify either "cp" or "mv".'
        exit(msg.format(action))

    if verbose:
        msg = 'Relocating {} to {}\n'
        stdout.write(msg.format(src_uri, dest_uri))
    subprocess.run(' '.join(['gsutil -m', action, src_uri, dest_uri]), shell=True)

