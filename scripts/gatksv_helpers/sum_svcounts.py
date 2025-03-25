#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Sum outputs from two or more runs of svtk count-svtypes
"""


import argparse
import csv
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('counts', help='One or more svtk count .tsv files', 
                        nargs='+')
    parser.add_argument('-o', '--outfile', help='output .tsv [default: stdout]',
                        default='stdout')
    args = parser.parse_args()

    # Iterate over each input count file and sum in a nested dict
    counter = {}
    for count_file in args.counts:
        with open(count_file) as fin:
            for sid, svtype, k in csv.reader(fin, delimiter='\t'):
                # Skip header
                if not k.isdigit():
                    continue
                # Otherwise, update counter dict
                if svtype not in counter.keys():
                    counter[svtype] = dict()
                if sid not in counter[svtype].keys():
                    counter[svtype][sid] = 0
                counter[svtype][sid] += int(k)

    # Open connection to output file
    if args.outfile in 'stdout /dev/stdout -'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Write summed counts to outfile
    outfile.write('sample\tsvtype\tcount\n')
    for svtype in sorted(counter.keys()):
        for sid in sorted(counter[svtype].keys()):
            outfile.write('{}\t{}\t{}\n'.format(sid, svtype, counter[svtype][sid]))

    # Clear buffer
    outfile.close()


if __name__ == '__main__':
    main()
