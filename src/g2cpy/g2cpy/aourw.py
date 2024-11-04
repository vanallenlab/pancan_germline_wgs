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


