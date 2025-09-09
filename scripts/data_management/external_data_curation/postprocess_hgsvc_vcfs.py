#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Clean up HGSVC WGS variant calls for G2C QC
"""


import argparse
import numpy as np
import pysam
import sys
from collections.abc import Iterable


def clean_gts(record):
    """
    Convert phased + unphased GTs to standard unphased biallelic notation
    """

    for s in record.samples:

        gts = list(record.samples[s].get('GT'))
        ploidy = len(gts) - 1
        
        pgts = [g for g in gts[:-1] if g is not None]
        ugt = gts[-1]
        
        if len(pgts) == 0:
            cgts = [ugt] * ploidy
        
        else:
            if len(pgts) < ploidy:
                pgts.append(ugt)
            n_alt = len([g for g in pgts if g is not None and g > 0])
            n_ref = len([g for g in pgts if g is not None and g == 0])
            n_na = np.max([0, ploidy - (n_alt + n_ref)])
            cgts = ([1] * n_alt) + ([0] * n_ref) + ([None] * n_na)
        
        new_gts = tuple(sorted(cgts[0:ploidy], key=lambda x: (x is not None, x)))
        record.samples[s]['GT'] = new_gts

    return record


def format_sv(record, minimal=False):
    """
    Format SVs to loose GATK-SV convention (always positive SVLEN and always specify END)
    This requires replacing the original record with a new record due to being 
    unable to set END for records read from a file with pysam
    """

    svtype = record.info.get('SVTYPE')
    if isinstance(svtype, Iterable) \
    and not isinstance(svtype, str):
        svtype = svtype[0]
    if svtype == 'UNK':
        svtype = 'CPX'
    svlen = record.info.get('SVLEN', 0)
    if isinstance(svlen, Iterable):
        svlen = int(svlen[0])

    record.info['SVLEN'] = int(np.abs(svlen))

    # Make new record as a near identical copy of the original
    if minimal:
        newrec = record
    else:
        if svtype == 'DEL':
            svlen = record.info['SVLEN']
        elif svtype == 'INS':
            svlen = 1
        newrec = record.copy()
        if len(newrec.ref) > 1:
            newrec.ref = 'N'
        if not newrec.alts[0].startswith('<'):
            newrec.alts = ('<{}>'.format(svtype), )
        newrec.stop = newrec.pos + svlen

    return newrec


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument('vcf_out', help='Output vcf. Also accepts "stdout" and "-".')
    parser.add_argument('--sv', action='store_true', help='Perform extra steps ' +
                        'for SV curation [default: only clean up GTs]')
    parser.add_argument('--minimal-sv', action='store_true', help='Only enforce ' +
                        'non-negative SVLEN values (i.e., do not update ' +
                        'variant coordinates)')
    parser.add_argument('--no-gt', action='store_true', help='Disable GT cleanup')
    args = parser.parse_args()

    # Open connections to input VCF
    if args.vcf_in in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf_in)

    # Prep output VCF
    if args.vcf_out in '- stdout'.split():
        fout = pysam.VariantFile(sys.stdout, 'w', header=vcf.header)
    else:
        fout = pysam.VariantFile(args.vcf_out, 'w', header=vcf.header)

    # Process each record from input VCF
    for rec in vcf:

        # Copy record to allow edits
        newrec = rec.copy()

        # Clean up genotypes
        if not args.no_gt:
            newrec = clean_gts(newrec)

        # Reformat SVs
        if args.sv or args.minimal_sv:
            newrec = format_sv(newrec, minimal=args.minimal_sv)

        # Write cleaned record to file
        fout.write(newrec)

    fout.close()


if __name__ == '__main__':
    main()


