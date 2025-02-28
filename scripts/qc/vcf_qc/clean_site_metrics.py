#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Clean site-level QC metrics
"""


import argparse
import csv
import gzip
import numpy as np
from Bio import bgzf
from g2cpy import classify_variant
from os import remove
from re import sub
from sys import stdin


# Declare common constants
# Note that short variant subclasses are lowercase and structural variant
# subclasses are uppercase to distinguish them by size
var_classes = 'snv indel sv'.split()
var_subclasses = {'snv' : ['snv'],
                  'indel' : 'ins del'.split(),
                  'sv' : 'DEL DUP CNV INS INV CPX CTX'.split()}


def make_size_bins(bins_per_log10, max_size):
    """
    Create bin boundaries for compressed size distribution counter
    
    Bin ranges are defined as the inclusive lower bound for each bin, which
    stops exclusively at the incrementally next larger bin lower bound

    Note we use integer counts for size ~ [0, 9], then default to
    --size-bins-per-log10 after that for a cleaner representation of indels
    """

    size_steps = np.arange(0, 1, 1 / bins_per_log10)
    size_log_ceil = int(np.ceil(np.log10(max_size)) + 1)
    size_ge = list(np.concatenate([10 ** (o + size_steps) for o in range(1, size_log_ceil)]))
    return np.array([int(round(x, 0)) for x in list(range(0, 10)) + size_ge if x <= max_size])


def make_af_bins(bins_per_log10, sample_size):
    """
    Create AF bin boundaries for compressed AF distribution counter

    Note that here each bin is defined in the opposite direction as size;
    namely, the count of variants with AF less than the threshold, down to
    the threshold defined by the previous bin
    
    Same as for size bins, we use integer counts for AC ~ [1, 9]
    then revert to log-uniformly spaced steps according to --bins-per-log10-af
    """

    max_an = 2 * sample_size
    af_steps = np.arange(0, -1, -1 / bins_per_log10)
    af_log_floor = int(np.floor(np.log10(1 / max_an)))
    af_lt = list(np.concatenate([10 ** (o + af_steps) for o in range(0, af_log_floor, -1)]))
    af_lt.reverse()
    smallest_nonac_bin = 10 / max_an
    ac_bins = [ac / max_an for ac in range(2, 10)]
    af_lt = ac_bins + [k for k in af_lt if k >= smallest_nonac_bin]
    af_lt[-1] = af_lt[-1] + 10e-6
    return np.array(af_lt)


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', help='Tab-separated QC metrics',
                        default='stdin')
    parser.add_argument('-o', '--output-prefix', help='Prefix for all output files',
                        metavar='[path]', default='./vcf_qc')
    parser.add_argument('-N', '--sample-count', type=int, required=True,
                        help='Effective number of samples in callset; required ' +
                        'for determining AF bins', metavar='[int]')
    parser.add_argument('--bins-per-log10-af', default=10, type=int,
                        help='Number of discrete bins to use per order of ' +
                        'magnitude when summarizing AF distributions',
                        metavar='[int]')
    parser.add_argument('--bins-per-log10-size', default=10, type=int,
                        help='Number of discrete bins to use per order of ' +
                        'magnitude when summarizing size distributions',
                        metavar='[int]')
    parser.add_argument('--max-size-bin', default=10e5, type=int,
                        help='Maximum size for binning; all variants larger than ' +
                        'this size will be included in the top bin.',
                       metavar='[int]')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress ' +
                        'output files with gzip/bgzip [default: write ' + 
                        'uncompressed .tsv/.bed]')
    args = parser.parse_args()

    # Create counter for compressed variant size distribution
    size_ge = make_size_bins(args.bins_per_log10_size, args.max_size_bin)
    size_counter = {vc : {vsc : [0] * len(size_ge) for vsc in var_subclasses[vc]}
                     for vc in var_classes}

    # Create counter for compressed AF distribution
    af_lt = make_af_bins(args.bins_per_log10_af, args.sample_count)
    af_counter = {vc : {vsc : [0] * len(af_lt) for vsc in var_subclasses[vc]}
                  for vc in var_classes}

    # Open connection to input file
    if args.input in '- stdin /dev/stdin'.split():
        fin = stdin
    else:
        fin = open(args.input)
    indat = csv.reader(fin, delimiter='\t')

    # Open connections to output BEDs with all site information
    site_fouts = {}
    site_header = '#chrom start end vid class subclass size ac af hwe exhet'.split()
    for vc in var_classes:
        fnbase = args.output_prefix + '.' + vc + '.sites.bed'
        if args.gzip:
            fout = bgzf.BgzfWriter(fnbase + '.gz', mode='wt')
        else:
            fout = open(fnbase, mode='w')
        fout.write('\t'.join(site_header) + '\n')
        site_fouts[vc] = fout

    # Iterate over each line in metrics and process variants one at a time
    vc_seen = {vc : False for vc in var_classes}
    for chrom, pos, end, ref, alt, svlen, ac, af, hwe, exhet in indat:

        # Clean up length based on ref/alt and svlen
        if svlen != '.':
            varlen = int(svlen)
        else:
            varlen = np.abs(len(alt) - len(ref))

        # Assign each variant to class & subclass
        vc, vsc = classify_variant(ref, alt, varlen)
        vc_seen[vc] = True

        # Assign simple VID
        if vc == 'sv':
            vid = '_'.join([str(x) for x in [chrom, pos, vsc, varlen]])
        else:
            vid = '_'.join([str(x) for x in [chrom, pos, ref, alt]])

        # Add variant count to binned size by class & subclass
        size_counter[vc][vsc][np.argmin(varlen >= size_ge)-1] += 1

        # Add variant count to binned AF by class & subclass
        af = float(af)
        af_counter[vc][vsc][np.argmax(af < af_lt)] += 1

        # Convert long floats to scientific notation
        af = '{:.2e}'.format(af)
        hwe = '{:.2e}'.format(float(hwe))
        exhet = '{:.2e}'.format(float(exhet))

        # Write reformatted variant data to out .bed
        outvals = [chrom, pos, end, vid, vc, vsc, varlen, ac, af, hwe, exhet]
        lout = '\t'.join([str(x) for x in outvals])
        site_fouts[vc].write(lout + '\n')

    # Close connections to all site-level output files to flush buffer
    for fc in site_fouts.values():
        fc.close()

    # Delete site-level summary files for which no qualifying variants were observred
    for vc, seen in vc_seen.items():
        if not seen:
            if args.gzip:
                delpath = site_fouts[vc]._handle.name
            else:
                delpath = site_fouts[vc].name
            remove(delpath)

    # Write compressed variant size distributions to file
    fnbase = args.output_prefix + '.size_distrib.tsv'
    if args.gzip:
        fout = gzip.open(fnbase + '.gz', mode='wt', encoding='utf-8')
    else:
        fout = open(fnbase, mode='w')
    header = '#class subclass'.split() + ['ge{}bp'.format(k) for k in size_ge]
    fout.write('\t'.join(header) + '\n')
    for vc, vsc_dat in size_counter.items():
        for vsc, counts in vsc_dat.items():
            if sum(counts) > 0:
                outvals = [vc, vsc] + [str(k) for k in counts]
                fout.write('\t'.join(outvals) + '\n')
    fout.close()

    # Write compressed AF distributions to file
    fnbase = args.output_prefix + '.af_distrib.tsv'
    if args.gzip:
        fout = gzip.open(fnbase + '.gz', mode='wt', encoding='utf-8')
    else:
        fout = open(fnbase, mode='w')
    header = '#class subclass'.split() + ['lt{:.2e}'.format(k) for k in af_lt]
    fout.write('\t'.join(header) + '\n')
    for vc, vsc_dat in af_counter.items():
        for vsc, counts in vsc_dat.items():
            if sum(counts) > 0:
                outvals = [vc, vsc] + [str(k) for k in counts]
                fout.write('\t'.join(outvals) + '\n')
    fout.close()


if __name__ == '__main__':
    main()

