#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Transform G2C VCF QC site statistics into approximate variant counts per sample
"""


# Import libraries
import argparse
import csv
import gzip
import numpy as np
from g2cpy import determine_filetype
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--sites-bed', default='stdin',
                        help='G2C QC-formatted sites.bed. Can be compressed. ' +
                        'Also accepts streamed standard input. [default: stdin]')
    parser.add_argument('--key-by', choices='vc vsc'.split(), default='vc',
                        help='Specify how counts should be aggregated. `vc` = ' +
                        'variant class, `vsc` = variant subclass [default: vc]')
    parser.add_argument('-o', '--outfile', default='stdout',
                        help='Path to output file [default: stdout]')
    args = parser.parse_args()

    # Open connection to input file
    if args.sites_bed in 'stdin /dev/stdin -'.split():
        handle = csv.reader(stdin, delimiter='\t')
    else:
        if 'compressed' in determine_filetype(args.sites_bed):
            fin = gzip.open(args.sites_bed, 'rt')
        else:
            fin = open(args.sites_bed)
        handle = csv.reader(fin, delimiter='\t')

    # Iterate over each line in --sites-bed
    res = {'all' : 0}
    for chrom, start, end, vid, vc, vsc, size, ac, af, fhet, fhom, hwe in handle:
        # Skip header line
        if chrom.startswith('#'):
            continue

        # Add key to result aggregator if not already present
        if args.key_by == 'vc':
            key = str(vc)
        elif args.key_by == 'vsc':
            key = str(vsc)
        if key not in res.keys():
            res[key] = 0

        # Estimate number of samples genotyped at this site
        ac = float(ac)
        af = float(af)
        an = ac / af
        nsamp = an / 2

        # Infer genotype frequencies according to HWE if not explicitly specified
        try:
            fhom = float(fhom)
        except:
            fhom = af ** 2
        try:
            fhet = float(fhet)
        except:
            fhet = 2 * af * (1 - af)

        # Estimate carrier rate for this variant
        ncarriers = nsamp * (fhom + fhet)
        vpg = ncarriers / nsamp
        res[key] += vpg
        res['all'] += vpg
    
    # Write results to output file
    if args.outfile in 'stdout /dev/stdout -'.split():
        fout = stdout
    else:
        fout = open(args.outfile, 'w')
    fout.write('#{}\tvariants_per_genome\n'.format(args.key_by))
    for k, v in sorted(res.items()):
        fout.write('{}\t{}\n'.format(k, int(np.round(v, 0))))
    fout.close()


if __name__ == '__main__':
    main()

