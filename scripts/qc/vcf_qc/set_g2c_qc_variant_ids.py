#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Name all variants in a VCF using G2C QC convention
"""


import argparse
import numpy as np
import pysam
import sys
from g2cpy import classify_variant, name_variant


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--vcf-in', default='stdin', help='Input vcf. Also ' +
                        'accepts "stdin" and "-". [default: stdin]')
    parser.add_argument('--vcf-out', default='stdout', help='Output vcf. Also ' +
                        'accepts "stdout" and "-". [default: stdout]')
    args = parser.parse_args()

    # Open connections to input VCF
    if args.vcf_in in '- stdin /dev/stdin'.split():
        invcf = pysam.VariantFile(sys.stdin)
    else:
        invcf = pysam.VariantFile(args.vcf_in)

    # Open connections to output VCF
    if args.vcf_out in '- stdout /dev/stdout'.split():
        outvcf = pysam.VariantFile(sys.stdout, mode='w', header=invcf.header)
    else:
        outvcf = pysam.VariantFile(args.vcf_out, mode='w', header=invcf.header)

    # Rename each record from --vcf-in before writing to --vcf-out
    for record in invcf:

        # Gather site descriptives
        ref, alt = record.alleles[:2]
        if 'SVLEN' in record.info.keys():
            varlen = int(record.info['SVLEN'])
        else:
            varlen = np.abs(len(alt) - len(ref))
        vc, vsc = classify_variant(ref, alt, varlen)
        vid = name_variant(record.chrom, record.pos, ref, alt, vc, vsc, varlen)
        record.id = vid

        # Write to --vcf-out
        outvcf.write(record)

    # Clear buffer
    invcf.close()
    outvcf.close()


if __name__ == '__main__':
    main()

