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


def chromsort(contigs):
    """
    Sorts an iterable of human contigs by their positional order
    """

    contig_order = {str(k + 1) : k + 1 for k in range(22)}
    contig_order.update({'X' : 23, 'Y' : 24})

    return sorted(contigs, key=lambda k: contig_order[sub('^chr', '', str(k))])


def classify_variant(ref, alt, var_len=None):
    """
    Maps a variant to its G2C variant class & subclass codes
    Returns two values: class, subclass
    """

    nucs = 'ACTG'.split()

    # Compute variant length if it is not already defined
    if var_len is None:
        var_len = abs(len(ref) - len(alt))

    # SNVs have zero size and cannot contain multinucleotide refs or alts
    if var_len == 0 and len(ref) == 1 and len(alt) == 1:
        return 'snv snv'.split()

    # Indels are 1-49bp
    elif var_len < 50:
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

    # SVs are >=50bp
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

