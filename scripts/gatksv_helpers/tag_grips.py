#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Tag deletions predicted to be involved in gene retroduplication events (gene retrocopy insertion polymorphism; GRIP)
"""


import argparse
import pybedtools as pbt
import pysam
import time
from sys import stdin, stdout


def make_intron_bed(gtf_in):
    """
    Build a map of all introns and return as pbt.BedTool
    """

    tx_info = {}
    prev_chrom = None
    intron_strs = ''
    ic_fmt = '{}\t{}\t{}\t{}\n'

    for record in pbt.BedTool(gtf_in):
        if record.fields[2] != 'exon':
            continue

        txid = record.attrs.get('transcript_id', None)
        if txid is None:
            continue

        # Check if we have moved to a new chromosome, in which case all old info
        # should be cleared from memory
        if prev_chrom is not None:
            if record.chrom != prev_chrom:
                tx_info = {}
        prev_chrom = record.chrom


        # If this transcript has been seen before, infer intron coordinates and add to intron_strs
        if txid in tx_info.keys():
            intron_coords = ic_fmt.format(record.chrom, tx_info[txid], 
                                          record.start, txid)
            intron_strs += intron_coords

        # Update latest exon junction
        tx_info[txid] = record.end

    return pbt.BedTool(intron_strs, from_string=True)


def intron_check(record, introns, ro=0.90):
    """
    Check if a record has reciprocal overlap with any introns
    """

    svc_fmt = '{}\t{}\t{}\n'
    rec_int = svc_fmt.format(record.chrom, record.start, record.stop)
    rec_bt = pbt.BedTool(rec_int, from_string=True)
    if len(introns.intersect(rec_bt, r=True, f=ro)) > 0:
        return True
    else:
        return False


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='input .vcf [default: stdin]')
    parser.add_argument('vcf_out', help='output .vcf [default: stdout]')
    parser.add_argument('--gtf', required=True, help='input .gtf')
    parser.add_argument('--pace', help='Report pace of progress', 
                        action='store_true')
    args = parser.parse_args()

    # Open connection to input vcf
    if args.vcf_in in 'stdin /dev/stdin -'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)

    # Add new FILTER tag to header
    header = invcf.header
    new_filt = '##FILTER=<ID=PREDICTED_GRIP_JXN,Description="This variant is ' + \
               'predicted to mark a splice junction for a gene retrocopy ' + \
               'insertion event and should not be evaluated as a canonical deletion.">'
    header.add_line(new_filt)

    # Build map of all introns
    introns = make_intron_bed(args.gtf)

    # Open connection to output VCF
    if args.vcf_out in 'stdout /dev/stdout -'.split():
        outvcf = pysam.VariantFile(stdout, mode='w', header=header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, mode='w', header=header)

    # Iterate over records in input VCF and mark deletions meeting GRIP criteria
    k = 0
    start = time.time()
    prev = start
    for record in invcf.fetch():
        if record.info['SVTYPE'] == 'DEL':
            if 'RD' not in record.info['EVIDENCE']:
                if intron_check(record, introns):
                    record.filter.add('PREDICTED_GRIP_JXN')
        outvcf.write(record)

        # Check pace if optioned
        k += 1
        if k % 10 == 0 and args.pace:
            now = time.time()
            print('Progress: checked {:,} variants in {:.1f} seconds (last 10 variants: {:.1f} secs)'.format(k, now - start, now - prev))
            prev = now

    # Clear buffer
    outvcf.close()


if __name__ == '__main__':
    main()
