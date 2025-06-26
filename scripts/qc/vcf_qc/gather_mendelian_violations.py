#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Collect Mendelian violation data from parent-child trios
"""


import argparse
import csv
import numpy as np
import pysam
import sys
from g2cpy import classify_variant, name_variant


mv_labels = 'denovo impossible_hom obligate_carrier obligate_hom pass'.split()


def load_trios(trios_in, all_samples):
    """
    Builds a dictionary of complete trios included in all_samples
    """

    trio_dict = {}

    with open(trios_in) as fin:
        for fid, proid, faid, moid, sex, pheno in csv.reader(fin, delimiter='\t'):
            if len(set([proid, faid, moid]).intersection(set(all_samples))) == 3:
                trio_dict[fid] = {'pro' : proid, 'fa' : faid, 'mo' : moid}
    
    return trio_dict


def classify_gt(gt):
    """
    Classifies a GT tuple by zygosity
    """

    alleles = set([a for a in gt if a is not None])

    if len(alleles) == 0:
        return None
    else:
        if 0 in alleles:
            return 'ref'
        else:
            return 'hom'


def mendelian_eval(record, pro, fa, mo):
    """
    Evaluates Mendelian consistency for trio genotypes
    """

    # Gather genotypes for each sample in the tro
    pro_gt = classify_gt(record.samples[pro]['GT'])
    fa_gt = classify_gt(record.samples[fa]['GT'])
    mo_gt = classify_gt(record.samples[mo]['GT'])

    # Skip sites with incomplete trio genotypes
    if any(gt is None for gt in [pro_gt, fa_gt, mo_gt]):
        return None

    # Count number of ref, het, and hom parents
    n_ref_par = sum(g == 'ref' for g in [fa_gt, mo_gt])
    n_het_par = sum(g == 'het' for g in [fa_gt, mo_gt])
    n_hom_par = 2 - (n_ref_par + n_het_par)

    # Work through decision tree of possible GT configurations
    if pro_gt == 'ref':
        if n_hom_par == 0:
            return 'pass'
        elif n_hom_par == 2:
            return 'obligate_hom'
        else:
            return 'obligate_carrier'

    elif n_ref_par == 2:
        return 'denovo'

    elif pro_gt == 'het':
        if n_par_hom == 2:
            return 'obligate_hom'
        else:
            return 'pass'

    elif pro_gt == 'hom':
        if n_ref_par > 0:
            return 'impossible_hom'
        else:
            return 'pass'

    else:
        print('Unexpected genotype configuration. Debug now:')
        import pdb; pdb.set_trace()


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('vcf_in', help='Input vcf. Also accepts "stdin" and "-".')
    parser.add_argument('--trios-fam', required=True, 
                        help='.fam file of trios to analyze')
    parser.add_argument('--common-af', type=float, default=0.01, 
                        help='AF cutoff for common variants. [default: 0.01]')
    parser.add_argument('--eligible-vids', help='File listing variant IDs to ' +
                        'be considered for analysis. [default: analyze all ' +
                        'variants in vcf_in]')
    parser.add_argument('--summary-out', default='stdout', help='Path to ' +
                        'output .tsv containing summary counts of Mendelian ' +
                        'analysis for each trio. Accepts "stdout" and "-". ' +
                        '[default: stdout]')
    args = parser.parse_args()

    # Open connections to input VCF
    if args.vcf_in in '- stdin'.split():
        vcf = pysam.VariantFile(sys.stdin)
    else:
        vcf = pysam.VariantFile(args.vcf_in)

    # Load trios as dict and subset to samples present in VCF
    samples_in_vcf = [sid for sid in vcf.header.samples]
    trio_dict = load_trios(args.trios_fam, samples_in_vcf)

    # Load list of eligible variant IDs, if optioned
    if args.eligible_vids is not None:
        with open(args.eligible_vids) as evid_in:
            elig_vids = set([vid.rstrip() for vid in evid_in.readlines()])
    else:
        elig_vids = None

    # Tally results in nested dict
    res_dict = {fid : {} for fid in trio_dict.keys()}

    # Process each record from input VCF
    for record in vcf:

        # Gather site descriptives
        ref, alt = record.alleles[:2]
        if 'SVLEN' in record.info.keys():
            varlen = int(record.info['SVLEN'])
        else:
            varlen = np.abs(len(alt) - len(ref))
        vc, vsc = classify_variant(ref, alt, varlen)
        vid = name_variant(record.chrom, record.pos, ref, alt, vc, vsc, varlen)

        # Skip variant if not included in --eligible-vids, if optioned
        if elig_vids is not None:
            if vid not in elig_vids:
                continue

        # Assign variant to AF bin
        if 'AF' in record.info.keys():
            freq = float(record.info['AF'][0])
        elif 'CN_NONREF_FREQ' in record.info.keys():
            freq = float(record.info['CN_NONREF_FREQ'])
        else:
            freq = 0
        if freq < args.common_af:
            freq_bin = 'lt{:.2e}'.format(args.common_af)
        else:
            freq_bin = 'ge{:.2e}'.format(args.common_af)

        # Process genotypes for each trio
        for fid, members in trio_dict.items():

            mv_label = mendelian_eval(record, *members.values())

            if mv_label is None:
                continue
            
            # Increment counter
            if vc not in res_dict[fid].keys():
                res_dict[fid][vc] = {}
            if vsc not in res_dict[fid][vc].keys():
                res_dict[fid][vc][vsc] = {}
            if freq_bin not in res_dict[fid][vc][vsc].keys():
                res_dict[fid][vc][vsc][freq_bin] = {mv : 0 for mv in mv_labels}
            res_dict[fid][vc][vsc][freq_bin][mv_label] += 1

    # Format results as compressed distribution and write to --summary-out
    if args.summary_out in '/dev/stdout stdout -'.split():
        sum_out = sys.stdout
    else:
        sum_out = open(args.summary_out, 'w')
    sum_out.write('\t'.join('#family_id class subclass freq_bin'.split() + mv_labels) + '\n')
    for fid, fdat in res_dict.items():
        for vc, vcdat in fdat.items():
            for vsc, vscdat in vcdat.items():
                for fbin, fdat in vscdat.items():
                    mv_k = [str(fdat.get(mv, 0)) for mv in mv_labels]
                    outline = '\t'.join([fid, vc, vsc, fbin] + mv_k)
                    sum_out.write(outline + '\n')
    sum_out.close()


if __name__ == '__main__':
    main()

