#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Genomics-related helper functions
"""


from re import sub


def chrom2int(contig):
    """
    Converts a human primary contig string to an integer for simpler sorting
    """

    contig_order = {'chr' + str(k + 1) : k + 1 for k in range(22)}
    contig_order.update({'chrX' : 23, 'chrY' : 24, 'chrM' : 25})

    return contig_order[contig]


def chromsort(contigs):
    """
    Sorts an iterable of human contigs by their positional order
    """

    return sorted(contigs, key=lambda k: chrom2int(str(k)))


def classify_variant(ref, alt, var_len=None):
    """
    Maps a variant to its G2C variant class & subclass codes
    Returns two values: class, subclass
    """

    pur = 'A G'.split()
    pyr = 'C T'.split()
    nucs = pur + pyr

    # Compute variant length if it is not already defined
    if var_len is None:
        var_len = abs(len(ref) - len(alt))

    # SNVs have zero size and cannot contain multinucleotide refs or alts
    if var_len == 0 and len(ref) == 1 and len(alt) == 1:
        if ref in pur and alt in pur:
            sc = 'ti'
        elif ref in pyr and alt in pyr:
            sc = 'ti'
        else:
            sc = 'tv'
        return 'snv', sc

    # Indels are 1-49bp, excluding certain SV types with indefinite sizes
    elif var_len < 50 and 'CTX' not in alt:
        if 'INS' in alt or 'DUP' in alt:
            return 'indel ins'.split()
        elif 'DEL' in alt:
            return 'indel del'.split()
        elif len(ref) > len(alt):
            return 'indel del'.split()
        elif len(ref) < len(alt):
            return 'indel ins'.split()
        elif len(ref) == len(alt):
            return 'indel cpx'.split()
        else:
            msg = 'g2cpy classify_variant unable to classify suspected indel ' + \
                  'with ref {}, alt {}, and length {:,}'
            exit(msg.format(ref, alt, var_len))

    # SVs are >=50bp, except for rare types with indefinite sizes
    else:
        if alt.startswith('<') and alt.endswith('>'):
            sc = sub('<|>', '', alt).split(':')[0]
            return 'sv', sc
        elif len(ref) > len(alt):
            return 'sv DEL'.split()
        elif len(ref) < len(alt):
            return 'sv INS'.split()
        else:
            msg = 'g2cpy classify_variant unable to classify suspected SV ' + \
                  'with ref {}, alt {}, and length {:,}'
            exit(msg.format(ref, alt, var_len))

