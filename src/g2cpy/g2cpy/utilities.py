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
