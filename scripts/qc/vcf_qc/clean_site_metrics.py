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
import subprocess
from Bio import bgzf
from g2cpy import classify_variant, name_variant
from os import remove
from re import sub
from sys import stdin


# Declare common constants
# Note that short variant subclasses are lowercase and structural variant
# subclasses are uppercase to distinguish them by size
var_classes = 'snv indel sv'.split()
var_subclasses = {'snv' : 'ti tv'.split(),
                  'indel' : 'ins del'.split(),
                  'sv' : 'DEL DUP CNV INS INV CPX CTX BND'.split()}


def make_size_bins(bins_per_log10, max_size):
    """
    Create bin boundaries for compressed size distribution counter
    
    Bin ranges are defined as the inclusive lower bound for each bin, which
    stops exclusively at the incrementally next larger bin lower bound

    Fractional bp are not permitted, so in practice there will only be seven bins
    from 1bp - 9bp (1, 2, 3, 4, 5, 6, 8)
    """

    size_steps = np.arange(0, 1, 1 / bins_per_log10)
    size_log_ceil = int(np.ceil(np.log10(max_size)) + 1)
    size_ge = np.concatenate([10 ** (o + size_steps) for o in range(0, size_log_ceil)])
    size_ge = np.insert(np.unique(np.round(size_ge)), 0, 0)
    return np.array([int(x) for x in size_ge if x <= max_size])


def make_af_bins(bins_per_log10, sample_size, min_af_bin):
    """
    Create AF bin boundaries for compressed AF distribution counter

    Note that here each bin is defined in the opposite direction as size;
    namely, the count of variants with AF less than the threshold, down to
    the threshold defined by the previous bin

    The first bin break will always be min_af_bin, if provided
    """

    max_an = 2 * sample_size
    af_steps = np.arange(0, -1, -1 / bins_per_log10)
    af_log_floor = int(np.floor(np.log10(1 / max_an)))
    af_lt = list(np.concatenate([10 ** (o + af_steps) for o in range(0, af_log_floor, -1)]))
    af_lt.reverse()
    if min_af_bin is None:
        min_af_bin = 2 / max_an
    af_lt = [min_af_bin] + [k for k in af_lt if k > min_af_bin]
    af_lt[-1] = af_lt[-1] + 10e-6
    return np.unique(np.array(af_lt))


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input', default='stdin', metavar='file',
                        help='Tab-separated QC metrics [default: stdin]')
    parser.add_argument('-o', '--output-prefix', help='Prefix for all output files',
                        metavar='[path]', default='./vcf_qc')
    parser.add_argument('-N', '--sample-count', type=int, required=True,
                        help='Effective number of samples in callset; required ' +
                        'for determining AF bins', metavar='[int]')
    parser.add_argument('--bins-per-log10-af', default=10, type=int,
                        help='Number of discrete bins to use per order of ' +
                        'magnitude when summarizing AF distributions',
                        metavar='[int]')
    parser.add_argument('--min-af-bin', type=float, help='Smallest AF cutoff ' +
                        'for second-smallest AF bin [default: doubletons; e.g., ' + 
                        '2 / (2 * --sample-count)]', metavar='[float]')
    parser.add_argument('--bins-per-log10-size', default=10, type=int,
                        help='Number of discrete bins to use per order of ' +
                        'magnitude when summarizing size distributions',
                        metavar='[int]')
    parser.add_argument('--max-size-bin', default=1e6, type=int,
                        help='Maximum size for binning; all variants larger than ' +
                        'this size will be included in the top bin.',
                       metavar='[int]')
    parser.add_argument('--common-af', type=float, help='AF cutoff for common ' +
                        'variants. If provided, will generate separate sites.bed ' +
                        'output files restricted to common variants. [default: ' +
                        'do not generate common variant-only .bed outputs]',
                        metavar='[float]')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress ' +
                        'output files with gzip/bgzip [default: write ' + 
                        'uncompressed .tsv/.bed]')
    args = parser.parse_args()

    # Create counter for compressed variant size distribution
    size_ge = make_size_bins(args.bins_per_log10_size, args.max_size_bin)
    size_counter = {vc : {vsc : [0] * len(size_ge) for vsc in var_subclasses[vc]}
                     for vc in var_classes}

    # Create counter for compressed AF distribution
    af_lt = make_af_bins(args.bins_per_log10_af, args.sample_count, args.min_af_bin)
    af_counter = {vc : {vsc : [0] * len(af_lt) for vsc in var_subclasses[vc]}
                  for vc in var_classes}

    # Create nested counters for compressed size vs. AF distribution
    x_size_ge = np.insert(10 ** np.arange(0, np.log10(args.max_size_bin) + 1), 0, 0)
    x_af_lt = 10 ** np.arange(np.ceil(np.log10(af_lt[0])), 1, 1)
    x_af_lt[-1] = x_af_lt[-1] + 1e-6
    x_counter = {vc : {vsc : {size : [0] * len(x_af_lt) for size in x_size_ge} 
                       for vsc in var_subclasses[vc]} for vc in var_classes}

    # Open connection to input file
    if args.input in '- stdin /dev/stdin'.split():
        fin = stdin
    else:
        fin = open(args.input)
    indat = csv.reader(fin, delimiter='\t')

    # Open connections to output BEDs with all site information
    site_fouts = {}
    common_fouts = {}
    site_header = '#chrom start end vid class subclass size ac af freq_het freq_hom hwe'.split()
    for vc in var_classes:
        fnbase = args.output_prefix + '.' + vc + '.sites'
        if args.gzip:
            fout = bgzf.BgzfWriter(fnbase + '.bed.gz', mode='wt')
        else:
            fout = open(fnbase + '.bed', mode='w')
        fout.write('\t'.join(site_header) + '\n')
        site_fouts[vc] = fout
        if args.common_af is not None:
            if args.gzip:
                cfout = bgzf.BgzfWriter(fnbase + '.common.bed.gz', mode='wt')
            else:
                cfout = open(fnbase + '.common.bed', mode='w')
            cfout.write('\t'.join(site_header) + '\n')
            common_fouts[vc] = cfout

    # Iterate over each line in metrics and process variants one at a time
    vc_seen = {vc : False for vc in var_classes}
    common_seen = {vc : False for vc in var_classes}
    for chrom, pos, end, ref, alt, svlen, an, ac, af, cnc, cnf, ac_het, ac_hom, ac_hemi, hwe in indat:

        # Clean up length based on ref/alt and svlen
        if svlen != '.':
            varlen = int(svlen)
        else:
            varlen = np.abs(len(alt) - len(ref))

        # Clean up AC/AF for CNVs
        if cnf != '.' and cnc != '.':
            ac = int(cnc)
            af = float(cnf)

        # Assign each variant to class & subclass
        vc, vsc = classify_variant(ref, alt, varlen)
        vc_seen[vc] = True

        # Assign simple VID
        vid = name_variant(chrom, pos, ref, alt, vc, vsc, varlen)

        # Set size to arbitrarily large value for gross interchromosomal 
        # translocations with undefined sizes
        if vsc == 'CTX':
            varlen = int(10e7)

        # Add variant count to binned size by class & subclass
        size_counter[vc][vsc][np.argmin(varlen >= size_ge)-1] += 1

        # Add variant count to binned AF by class & subclass
        if af != '.':
            af = float(af)
            af_counter[vc][vsc][np.argmax(af < af_lt)] += 1
            if args.common_af is not None \
            and af >= args.common_af:
                common_seen[vc] = True
            
            # Add variant count to binned 2D size X AF by class & subclass
            x_counter[vc][vsc][x_size_ge[np.argmin(varlen >= x_size_ge)-1]][np.argmax(af < x_af_lt)] += 1

        # Infer genotype frequencies for diploid loci
        freq_het = 'NA'
        freq_hom = 'NA'
        if an != '.' and ac_het != '.' and ac_hom != '.' and ac_hemi != '.':
            an = int(an)
            ac_het = int(ac_het)
            ac_hom = int(ac_hom)
            ac_hemi = int(ac_hemi)
            an_bi = an - ac_hemi
            n_bi = an_bi / 2
            if n_bi > 0:
                freq_het = '{:.2e}'.format(ac_het / n_bi)
                freq_hom = '{:.2e}'.format((ac_hom / 2) / n_bi)

        # Convert long floats to scientific notation
        if af == '.':
            af = 'NA'
        else:
            af = '{:.2e}'.format(af)
        if hwe == '.':
            hwe = 'NA'
        else:
            hwe = '{:.2e}'.format(float(hwe))

        # Write reformatted variant data to out .bed
        outvals = [chrom, pos, end, vid, vc, vsc, varlen, ac, af, freq_het, freq_hom, hwe]
        lout = '\t'.join([str(x) for x in outvals])
        site_fouts[vc].write(lout + '\n')
        if args.common_af is not None \
        and af != '.':
            if float(af) >= args.common_af:
                common_fouts[vc].write(lout + '\n')

    # Close connections to all site-level output files to flush buffer
    for fc in list(site_fouts.values()) + list(common_fouts.values()):
        fc.close()

    # Delete site-level summary files for which no qualifying variants were observed
    for vc, seen in vc_seen.items():
        if args.gzip:
            fpath = site_fouts[vc]._handle.name
        else:
            fpath = site_fouts[vc].name
        
        # Delete file if empty
        if not seen:
            remove(fpath)

        # Otherwise, index summary files with tabix if --gzip is optioned
        else:
            if args.gzip:
                subprocess.run(['tabix', '-f', fpath], check=True)

    # Delete site-level common variant summary files for which no qualifying variants were observed
    if args.common_af is not None:
        for vc, seen in common_seen.items():
            if args.gzip:
                fpath = common_fouts[vc]._handle.name
            else:
                fpath = common_fouts[vc].name
            
            # Delete file if empty
            if not seen:
                remove(fpath)

            # Otherwise, index summary files with tabix if --gzip is optioned
            else:
                if args.gzip:
                    subprocess.run(['tabix', '-f', fpath], check=True)

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

    # Write compressed 2D size X AF distributions to file
    fnbase = args.output_prefix + '.size_vs_af_distrib.tsv'
    if args.gzip:
        fout = gzip.open(fnbase + '.gz', mode='wt', encoding='utf-8')
    else:
        fout = open(fnbase, mode='w')
    header = '#class subclass size'.split() + ['lt{:.2e}'.format(k) for k in x_af_lt]
    fout.write('\t'.join(header) + '\n')
    for vc, vsc_dat in x_counter.items():
        for vsc, size_dat in vsc_dat.items():
            for size, counts in size_dat.items():
                if sum(counts) > 0:
                    outvals = [vc, vsc, 'ge{:d}'.format(int(size))] + [str(k) for k in counts]
                    fout.write('\t'.join(outvals) + '\n')
    fout.close()


if __name__ == '__main__':
    main()

