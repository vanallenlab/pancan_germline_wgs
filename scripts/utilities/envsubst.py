#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Python implementation of bash envsubst
"""


# Import libraries
import argparse
from os.path import expandvars
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--infile', help='input file [default: stdin]', 
                        default='stdin')
    parser.add_argument('-o', '--outfile', help='output file [default: stdout]', 
                        default='stdout')
    args = parser.parse_args()

    # Open connection to input 
    if args.infile in '- stdin /dev/stdin'.split():
        fin = stdin
    else:
        fin = open(args.infile, 'r')

    # Open connection to output
    if args.outfile in '- stdout /dev/stdout'.split():
        fout = stdout
    else:
        fout = open(args.outfile, 'w')

    # Read lines one at a time from infile, substitute variables, and write to outfile
    for line in fin.readlines():
        fout.write(expandvars(line))

    # Close output file to clear buffer
    fout.close()


if __name__ == '__main__':
    main()
