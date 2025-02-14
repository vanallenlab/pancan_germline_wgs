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
