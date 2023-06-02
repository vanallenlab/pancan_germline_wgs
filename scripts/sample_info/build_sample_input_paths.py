#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Genomic Architecture of Human Cancers
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Build table of gs:// paths to expected inputs for samples from a single cohort in DFCI-G2C
"""


import argparse
import csv
from sys import stdout


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('cohort', help='name of cohort')
    parser.add_argument('samples', help='list of samples in cohort')
    parser.add_argument('-o', '--outfile', default='stdout',
                        help='output .tsv [default: stdout]')
    parser.add_argument('-b', '--bucket', default='gs://dfci-g2c-inputs',
                        help='Root bucket')
    parser.add_argument('--no-header', action='store_true', help='Do not ' +
                        'write header to --outfile')
    args = parser.parse_args()

    # Prepare expected path formatting conventions
    base_path = '{}/{}'.format(args.bucket, args.cohort)
    gs_fmts = {'gatksv_coverage' : '{}/gatk-sv/coverage/{}.counts.tsv.gz',
               'gatksv_pe' : '{}/gatk-sv/pesr/{}.pe.txt.gz',
               'gatksv_sd' : '{}/gatk-sv/pesr/{}.sd.txt.gz',
               'gatksv_sr' : '{}/gatk-sv/pesr/{}.sr.txt.gz',
               'gatkhc_gvcf' : '{}/gatk-hc/{}.g.vcf.gz',
               'manta_vcf' : '{}/manta/{}.manta.vcf.gz',
               'melt_vcf' : '{}/melt/{}.melt.vcf.gz',
               'wham_vcf' : '{}/wham/{}.wham.vcf.gz'}

    # Open connection to outfile
    if args.outfile in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.outfile, 'w')

    # Write header
    if not args.no_header:
        header = '\t'.join(['#sample_id'] + list(gs_fmts.keys())) + '\n'
        outfile.write(header)

    with open(args.samples) as fin:
        for line in fin.readlines():
            sid = line.rstrip()
            outvals = [sid] + [k.format(base_path, sid) for k in gs_fmts.values()]
            outfile.write('\t'.join(outvals) + '\n')

    # Close connection to otufile to clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

