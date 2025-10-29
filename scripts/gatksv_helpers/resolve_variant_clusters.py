#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Extract and collapse clusters of records from an input VCF
"""


import argparse
import csv
import numpy as np
import pysam
from g2cpy import recursive_flatten, integrate_infos, integrate_gts
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--in-vcf', metavar='vcf', default='stdin',
                        help='input .vcf [default: stdin]')
    parser.add_argument('--clusters', metavar='bed', help='BED4 with one line ' +
                        'per cluster to process. Fourth column is comma-delimited ' +
                        'list of cluster members', required=True)
    parser.add_argument('--out-vcf', metavar='vcf', default='stdout',
                        help='output .vcf [default: stdout]')
    parser.add_argument('-p', '--prefix', help='Name prefix for clustered variants',
                        default='reclustered')
    parser.add_argument('-b', '--buffer', type=int, default=10, 
                        help='Buffer for VCF queries, in base pairs')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.in_vcf in 'stdin /dev/stdin -'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.in_vcf)

    # Open connection to cluster BED file
    cluster_bed = csv.reader(open(args.clusters), delimiter='\t')

    # Get list of samples from input vcf
    samples = [s for s in invcf.header.samples]

    # Open connection to output vcf
    if args.out_vcf in 'stdout /dev/stdout -'.split():
        outvcf = pysam.VariantFile(stdout, mode='w', header=invcf.header)
    else:
        outvcf = pysam.VariantFile(args.out_vcf, mode='w', header=invcf.header)

    # Process each cluster in serial
    k = 0
    for chrom, start, end, mem_str in cluster_bed:
        
        vids = mem_str.split(',')

        # Gather records from VCF
        records = []
        for rec in invcf.fetch(chrom, 
                               np.nanmax([int(start) - args.buffer, 0]), 
                               int(end) + args.buffer):
            if rec.id in vids:
                records.append(rec)

        # Check that all records can be found
        if len(records) != len(vids):
            msg = 'Unable to find all {:,} records in --in-vcf for cluster at {}:{}-{}'
            exit(msg.format(len(vids), chrom, start, end))
        k += 1

        # Use first record as a template
        newrec = records[0].copy()

        # Assign basic (non-INFO) record information
        newrec.pos = int(np.floor(np.nanmedian([r.pos for r in records])))
        newrec.stop = int(np.floor(np.nanmedian([r.stop for r in records])))
        newrec.info['SVLEN'] = int(np.nanmax([newrec.stop - newrec.pos, 0]))
        newrec.id = '{}_{}'.format(args.prefix, k)
        newrec.qual = int(np.round(np.nanmean([r.qual for r in records])))
        newrec.filter.clear()
        for f in list(set(recursive_flatten([r.filter.items() for r in records]))):
            newrec.filter.add(f)

        # Merge INFOs
        newrec.info.clear()
        newrec.info.update(integrate_infos(records))

        # Merge GTs
        newrec = integrate_gts(newrec, records)

        # Write clustered record to output VCF
        outvcf.write(newrec)

    # Clear buffers
    outvcf.close()


if __name__ == '__main__':
    main()

