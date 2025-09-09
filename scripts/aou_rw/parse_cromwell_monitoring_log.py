#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Extract resource information from a single Cromwell monitoring.log
"""


# Import libraries
import argparse
import pandas as pd
import re
from os import path
from sys import stdout


def infer_task_name(logfile):
    """
    Infer task name/family from the path of the logfile
    """

    dn = path.dirname(logfile)
    parts = [p for p in dn.split('/') if p.startswith('call-')]
    return ':'.join([re.sub('^call-', '', p) for p in parts])



def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('log', help='Cromwell resource monitoring.log to parse')
    parser.add_argument('-o', '--outfile', default='stdout',
                        help='Path to output file [default: stdout]')
    parser.add_argument('--header', default=False, action='store_true',
                        help='Add header to --outfile')
    args = parser.parse_args()

    # Infer task family based on log name
    tname = infer_task_name(args.log)

    # Task resource allocations are printed in the header of the file
    # We need to scrape these first
    res_alloc = {'runtime' : 'NA'}
    with open(args.log) as fin:
        for line in fin:
            if line.startswith('Num processors:'):
                res_alloc['cpu'] = int(line.rstrip().split(': ')[1])
            if line.startswith('Total Memory:'):
                res_alloc['mem'] = float(line.rstrip().split(': ')[1].split()[0])
            if line.startswith('Total Disk space:'):
                res_alloc['disk'] = float(line.rstrip().split(': ')[1].split()[0])
            if 'Runtime Information' in line:
                header = fin.readline().rstrip().split()
                break

    # Read the body of the monitor log as a pd.DataFrame
    df = pd.read_csv(args.log, skiprows=8, sep='\t', header=None, names=header)

    # Compute peak stats
    h, m, s = map(int, df.ElapsedTime.iloc[-1].split(':'))
    peak = {
        'mem' : df.Mem.max(),
        'disk' : df.Disk.max(),
        'cpu' : res_alloc.get('cpu', 1) * df.CPU.max() / 100,
        'runtime' : max([10, (3600 * h) + (60 * m) + s])
    }

    # Open connection to output file and optionally write header
    if args.outfile in 'stdout /dev/stdout -'.split():
        fout = stdout
    else:
        fout = open(args.outfile, 'w')
    if args.header:
        fout.write('\t'.join('#task resource allocated peak_used'.split()) + '\n')

    # Report results to output file
    for key, val in peak.items():
        alloc = res_alloc[key]
        outline = '\t'.join(map(str, [tname, key, alloc, val])) + '\n'
        fout.write(outline)
    fout.close()


if __name__ == '__main__':
    main()
