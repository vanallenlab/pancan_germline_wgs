#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Divide a GATK-style interval_list file into shards of a desired average size
"""


# Import libraries
import argparse
import numpy as np
import pandas as pd
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-intervals', required=True, 
                        help='Path to input .interval_list', default='stdin')
    parser.add_argument('-t', '--target-size', help='Desired average interval size',
                        type=float, default=10e10, required=True)
    parser.add_argument('-o', '--output-intervals', required=True,
                        help='Path to output .interval_list', default='stdin')
    args = parser.parse_args()

    # Open connection to input file
    if args.input_intervals in '- stdin /dev/stdin'.split():
        infile = stdin
    else:
        infile = open(args.input_intervals)

    # Open connection to output file
    if args.output_intervals in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.output_intervals, 'w')

    # Read input file into memory
    int_list = []
    for line in infile.readlines():

        # Write header from input to output
        if line.startswith('@'):
            outfile.write(line)

        # Store input intervals as tuples of (interval size, list of fields)
        else:
            fields = line.rstrip().split()
            size = int(fields[2]) - int(fields[1])
            int_list.append((size, fields, ))

    # Convert input data to two-column pd.DataFrame
    int_df = pd.DataFrame(int_list, columns='isize fields'.split())

    # Split the largest interval in half until mean interval is ≤ --target-size
    avg_size = int_df.isize.mean()
    while avg_size > args.target_size:

        # Find largest interval
        big_idx = int_df.isize.argmax()
        parent = int_df.fields[big_idx]
        parent[1] = int(parent[1])
        parent[2] = int(parent[2])

        # Divide interval in half
        midpoint = int(np.floor(np.mean([parent[2], parent[1]])))
        children = [(midpoint - parent[1], [*parent[0:2], midpoint, *parent[3:]], ),
                    (parent[2] - midpoint, [parent[0], midpoint+1, *parent[2:]], ),]

        # Drop parent and append children
        int_df = pd.concat([int_df.drop(big_idx),
                            pd.DataFrame(children, columns=int_df.columns)]).\
                    reset_index(drop=True)

        # Recompute mean interval size
        avg_size = int_df.isize.mean()

    # Sort sharded intervals and write to output file
    for vals in sorted(int_df.fields.tolist(), key=lambda x: int(x[1])):
        outfile.write('\t'.join([str(x) for x in vals]) + '\n')

    # Clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

