#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Extract peak LD values per variant class
"""


import argparse
import csv
import gzip
import numpy as np
from g2cpy import determine_filetype
from sys import stdout


def load_vid_list(infile):
    """
    Load a list of eligible variant IDs from a file and store as keys in a dict
    """

    if infile is not None:
        with open(infile) as fin:
            return {l.rstrip(): {} for l in fin.readlines()}
    else:
        return {}


def inject_nulls(res, null_r2):
    """
    Initialize each variant peak LD as a null value
    """

    has_snvs = len(res['snv']) > 0
    has_indels = len(res['indel']) > 0
    has_svs = len(res['sv']) > 0

    for vc in 'snv indel sv'.split():
        for vid in res[vc]:
            nr2 = null_r2.get(vid, 0)
            if has_snvs:
                res[vc][vid]['snv'] = nr2
            if has_indels:
                res[vc][vid]['indel'] = nr2
            if has_svs:
                res[vc][vid]['sv'] = nr2

    return res


def determine_vc(vid, snv={}, indel={}, sv={}):
    """
    Determines the class of a variant based on its ID
    """

    if vid in snv.keys():
        return 'snv'
    elif vid in indel.keys():
        return 'indel'
    elif vid in sv.keys():
        return 'sv'
    else:
        return None


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcor_slim', metavar='.tsv', 
                        help='Three-column .tsv of VID 1, VID 2, LD R2')
    parser.add_argument('--snv-list', help='List of eligible SNV variant IDs')
    parser.add_argument('--indel-list', help='List of eligible indel variant IDs')
    parser.add_argument('--sv-list', help='List of eligible SV variant IDs')
    parser.add_argument('--null-r2', metavar='.tsv', 
                        help='Two-column .tsv specifying custom null R2 values '+
                        'to use for each variant if no match is found.')
    parser.add_argument('-o', '--out-tsv', metavar='.tsv', 
                        help='Path to output .tsv [default: stdout]')
    args = parser.parse_args()

    # Check that at least one of SNV, indel, or SV IDs were provided
    if all([x is None for x in [args.snv_list, args.indel_list, args.sv_list]]):
        stop('At least one of --snv-list, --indel-list, or --sv-list must be provided.')

    # Load variant IDs and structure them as a results dictionary
    res = {
        'snv' : load_vid_list(args.snv_list),
        'indel' : load_vid_list(args.indel_list),
        'sv' : load_vid_list(args.sv_list)
    }

    # Load null R2 values, or fill with zeroes
    if args.null_r2 is not None:
        with open(args.null_r2) as fin:
            null_r2 = {vid : float(r2) for vid, r2 in csv.reader(fin, delimiter='\t')}
    else:
        null_r2 = {}
    res = inject_nulls(res, null_r2)

    # Iterate over each record in vcor-slim
    with open(args.vcor_slim) as fin:
        for v1, v2, r2 in csv.reader(fin, delimiter='\t'):
            vc1 = determine_vc(v1, *res.values())
            vc2 = determine_vc(v2, *res.values())
            if any([vc is None for vc in [vc1, vc2]]):
                continue
            res[vc1][v1][vc2] = np.nanmax([res[vc1][v1].get(vc2, 0), float(r2)])
            res[vc2][v2][vc1] = np.nanmax([res[vc2][v2].get(vc1, 0), float(r2)])

    # Open connection to output file
    if args.out_tsv in 'stdout /dev/stdout -'.split():
        fout = stdout
    else:
        if 'compressed' in determine_filetype(args.out_tsv):
            fout = gzip.open(args.out_tsv, 'wt')
        else:
            fout = open(args.out_tsv, 'w')

    # Write peak LD results to output file
    fout.write('\t'.join('#vid other_vc ld_r2'.split()) + '\n')
    for vc in 'snv indel sv'.split():
        for vid, subres in res[vc].items():
            for vc2, r2 in subres.items():
                fout.write('{}\t{}\t{:.4f}\n'.format(vid, vc2, r2))
    fout.close()


if __name__ == '__main__':
    main()

