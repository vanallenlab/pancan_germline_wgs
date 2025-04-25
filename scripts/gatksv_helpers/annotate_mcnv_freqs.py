#!/usr/bin/env python
# -*- coding: utf-8 -*-
#

"""
Compute copy-number frequencies for multiallelic CNVs

This is a greatly simplified version of GATK-SV compute_AFs.py
"""


import sys
import argparse
import pysam
from svtk import utils as svu
from collections import Counter


def calc_cn_freq(record, samples, ref_CN=2):
    """
    Computes allele frequencies for a single record based on a list of samples
    """

    # Get all sample CNs and remove Nones
    CNs_wNones = [record.samples[s]['CN'] for s in samples]
    CNs = [c for c in CNs_wNones if c is not None and c not in '. NA'.split()]

    if len(CNs) == 0:
        nonnull_CNs, nonref_CN_count, nonref_CN_freq = [0] * 3
        CN_dist = (0, )
        CN_freqs = (0, )
        CN_status = (0, )
    else:
        # Count number of samples per CN and total CNs observed
        CN_counts = dict(Counter(CNs))
        nonnull_CNs = len(CNs)

        # Get max observed CN and enumerate counts/frequencies per CN as list starting from CN=0
        max_CN = max([int(k) for k, v in CN_counts.items()])
        CN_dist = [int(CN_counts.get(k, 0)) for k in range(max_CN + 1)]
        CN_freqs = [round(v / nonnull_CNs, 6) for v in CN_dist]
        CN_status = [s for s in range(max_CN + 1)]

        # Get total non-reference CN counts and freq
        nonref_CN_count = sum([int(CN_counts.get(k, 0)) for k in range(max_CN + 1) if k != ref_CN])
        nonref_CN_freq = round(nonref_CN_count / nonnull_CNs, 6)

    # Add values to INFO field
    record.info['CN_NUMBER'] = nonnull_CNs
    record.info['CN_COUNT'] = tuple(CN_dist)
    record.info['CN_FREQ'] = tuple(CN_freqs)
    record.info['CN_STATUS'] = tuple(CN_status)
    record.info['CN_NONREF_COUNT'] = nonref_CN_count
    record.info['CN_NONREF_FREQ'] = nonref_CN_freq

    return record


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf', help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument('fout', help='Output vcf. Also accepts "stdout" and "-".')
    args = parser.parse_args()

    # Open connections to input VCF
    if args.vcf in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf)

    # Get list of all samples in vcf
    samples_list = list(vcf.header.samples)

    # Add relevant fields to header
    INFO_ADD = [
        '##INFO=<ID=CN_NUMBER,Number=1,Type=Integer,Description="Total number of samples with estimated copy numbers (multiallelic CNVs only).">',
        '##INFO=<ID=CN_COUNT,Number=.,Type=Integer,Description="Number of samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">',
        '##INFO=<ID=CN_STATUS,Number=.,Type=Integer,Description="Copy states corresponding to CN_COUNT, CN_FREQ: 0,1,...,maximum observed copy state (multiallelic CNVs only).">',
        '##INFO=<ID=CN_FREQ,Number=.,Type=Float,Description="Frequency of samples observed at each copy state, starting from CN=0 (multiallelic CNVs only).">',
        '##INFO=<ID=CN_NONREF_COUNT,Number=1,Type=Integer,Description="Number of samples with non-reference copy states (multiallelic CNVs only).">',
        '##INFO=<ID=CN_NONREF_FREQ,Number=1,Type=Float,Description="Frequency of samples with non-reference copy states (multiallelic CNVs only).">'
    ]

    for line in INFO_ADD:
        vcf.header.add_line(line)

    # Prep output VCF
    if args.fout in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.fout, 'w', header=vcf.header)

    # Iterate over VCF and annotate multiallelic records
    for r in vcf:
        if not svu.is_biallelic(r):
            r = calc_cn_freq(r, samples_list)
        fout.write(r)

    fout.close()


if __name__ == '__main__':
    main()

