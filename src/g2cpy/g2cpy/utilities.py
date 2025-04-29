#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Generic small / utility functions
"""


def astype_default(x, type_fxn, default=None):
    """
    Performs "soft" fault-tolerant type coercion of `x` to type `type_fxn()`
    Any values that fail to be coerced by `type_fxn()` will be assigned `default`
    """

    try:
        return type_fxn(x)
    except (ValueError, TypeError):
        return type_fxn(default)


def recursive_flatten(S):
    """
    Recursively flatten all elements of an input list, S
    Taken from SO: https://stackoverflow.com/questions/12472338/flattening-a-list-recursively
    """

    if S == []:
        return S
    if isinstance(S[0], list):
        return recursive_flatten(S[0]) + recursive_flatten(S[1:])
    return S[:1] + recursive_flatten(S[1:])
