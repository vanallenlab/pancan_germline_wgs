#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Bifurcate an SV VCF to prepare for SNV-based GT refinement
"""


import argparse
import pysam


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-vcf', required=True, metavar='VCF', 
                        type=str, help='Input .vcf. Required.')
    parser.add_argument('-e', '--eligible-output-vcf', 
                        help='Output .vcf for eligible SVs',
                        required=True, metavar='VCF', type=str)
    parser.add_argument('-p', '--passthrough-output-vcf', 
                        help='Output .vcf for ineligible "pass-through" SVs',
                        metavar='VCF', type=str)
    parser.add_argument('--min-af', default=0.05, type=float, metavar='Float',
                        help='Lower AF threshold for eligibility [default: 0.05]')
    parser.add_argument('--max-af', default=0.95, type=float, metavar='Float',
                        help='Upper AF threshold for eligibility [default: 0.95]')
    parser.add_argument('--min-ac', default=20, type=int, metavar='Int',
                        help='Minimum AC for eligibility [default: 20]')
    parser.add_argument('--max-ac', default=10e10, type=int, metavar='Int',
                        help='Minimum AC for eligibility [default: ~infinite]')
    parser.add_argument('--samples', type=chr, metavar='.txt',
                        help='List of sample IDs to include [default: keep all samples]')
    args = parser.parse_args()

    # Open connection to input VCF
    invcf = pysam.VariantFile(args.input_vcf)
    all_samples = set([s for s in invcf.header.samples])

    # Modify output header to only include --samples, if optioned
    # TODO: implement this

    # Open connection to output VCF(s)
    # TODO: implement this

    # Iterate over input VCF and route each SV to the correct output VCF
    # TODO: implement this

    # Close connections to output VCF(s)
    # TODO: implement this

if __name__ == '__main__':
    main()

