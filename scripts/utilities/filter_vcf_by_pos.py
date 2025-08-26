#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Extract records from a VCF based on strict POS-based filtering
This differs from default interval overlap-based filtering for htslib/bedtools
Only records with POS contained within --regions bed will be included in --out-vcf
"""


# Import libraries
import argparse
import pandas as pd
import pybedtools as pbt
import pysam
from g2cpy import name_record
from os import path
from sys import stdin, stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-vcf', default='stdin',
                        help='Input VCF. Index-jumps instead of streaming if ' +
                        '.tbi is able to be found [default: stdin]')
    parser.add_argument('-r', '--regions-bed', required=True,
                        help='BED file of regions to extract')
    parser.add_argument('-g', '--genome-file', help='BEDTools-style .genome ' +
                        'file. If provided, will be used to sort --regions-bed, ' +
                        'which will otherwise be sorted with BEDTools default ' +
                        'ordering')
    parser.add_argument('-o', '--out-vcf', default='stdout',
                        help='Output VCF [default: stdout]')
    args = parser.parse_args()

    # Load target regions as pd.DataFrame
    tbt = pbt.BedTool(args.regions_bed).cut(range(3))
    if args.genome_file is not None:
        tbt = tbt.sort(g=args.genome_file).merge()
    else:
        tbt = tbt.sort().merge()
    tdf = tbt.to_dataframe()

    # Open connection to input VCF
    ijump = False
    if args.input_vcf in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(stdin)
    else:
        invcf = pysam.VariantFile(args.input_vcf)
        if path.exists(args.input_vcf + '.tbi'):
            ijump = True

    # Open connection to output VCF
    if args.out_vcf in '- stdout /dev/stdout'.split():
        outvcf = pysam.VariantFile(stdout, 'w', header=invcf.header)
    else:
        outvcf = pysam.VariantFile(args.out_vcf, 'w', header=invcf.header)

    # Extract variants based on value of ijump
    vids_seen = set()
    if ijump:
        for region in tbt:
            for rec in invcf.fetch(region.chrom, region.start, region.end):
                if rec.id == '.' or rec.id is None:
                    rec.id = name_record(rec)
                if rec.pos >= region.start \
                and rec.stop <= region.end \
                and rec.id not in vids_seen:
                    outvcf.write(rec)
                    vids_seen.add(rec.id)
    else:

    # Clear buffer unless writing to stdout
    if outvcf != stdout:
        outvcf.close()


if __name__ == '__main__':
    main()

