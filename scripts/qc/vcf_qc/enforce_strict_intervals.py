#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Enforce strict adherence to query intervals for a BED file of variant summary metrics
"""


# Import libraries
import argparse
import gzip
import pandas as pd
import pybedtools as pbt
from g2cpy import determine_filetype


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-bed', required=True, 
                        help='BED file of variants to filter')
    parser.add_argument('-t', '--target-intervals', required=True,
                        help='BED file of intervals for filtering')
    parser.add_argument('-f', '--frac-coverage', type=float, default=0.5,
                        help='Minimum fraction of each variant that must be ' +
                             'covered by --target-intervals')
    parser.add_argument('-m', '--min-size', type=int, default=0,
                        help='Minimum variant size')
    parser.add_argument('-M', '--max-size', type=int, default=10e10,
                        help='Maximum variant size')
    parser.add_argument('-o', '--output-bed', required=True,
                        help='Path to output BED of filtered variants')
    args = parser.parse_args()

    # Load target intervals as both pbt.BedTool and pd.DataFrame
    tbt = pbt.BedTool(args.target_intervals)
    tdf = tbt.to_dataframe()

    # Read and save header line from --input-bed
    if 'compressed' in determine_filetype(args.input_bed):
        with gzip.open(args.input_bed,  'rt') as fin:
            header = fin.readline().rstrip()
    else:
        with open(args.input_bed) as fin:
            header = fin.readline().rstrip()
    if not header.startswith('#'):
        header = None

    # Open connection to input variants
    vbt = pbt.BedTool(args.input_bed)
    n_fields = len(vbt[0].fields)

    # Apply filters to variants and save to --output-bed
    vbt.filter(lambda x: int(x[6]) >= args.min_size and int(x[6]) <= args.max_size).\
        filter(lambda x: any((x.start >= tdf.start) & (x.start <= tdf.end))).\
        coverage(b=tbt).\
        filter(lambda x: float(x.fields[-1]) >= args.frac_coverage).\
        cut(range(n_fields)).\
        saveas(args.output_bed, trackline=header)


if __name__ == '__main__':
    main()


