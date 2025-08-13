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


import hashlib
import math


base62_alphabet = '0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ'


def astype_default(x, type_fxn, default=None):
    """
    Performs "soft" fault-tolerant type coercion of `x` to type `type_fxn()`
    Any values that fail to be coerced by `type_fxn()` will be assigned `default`
    """

    try:
        return type_fxn(x)
    except (ValueError, TypeError):
        return type_fxn(default)


def baseX_encode(num, alphabet=base62_alphabet):
    """
    Encodes an integer in a custom Base alphabet
    Taken partially from https://stackoverflow.com/questions/1119722/base-62-conversion
    """
    
    if num == 0:
        return alphabet[0]

    arr = []
    arr_append = arr.append
    _divmod = divmod
    base = len(alphabet)

    while num:
        num, rem = _divmod(num, base)
        arr_append(alphabet[rem])
    
    arr.reverse()
    
    return ''.join(arr)


def hash_string(input, out_length, alphabet=base62_alphabet):
    """
    Custom encoding & hashing function for any input string
    Designed to generate variant ID suffixes, but can theoretically be used for any string
    """

    bits_needed = math.ceil(out_length * math.log2(len(alphabet)))
    bytes_needed = math.ceil(bits_needed / 8)

    digest =  hashlib.sha256(input.encode('utf-8')).digest()[:bytes_needed]
    num = int.from_bytes(digest, 'big')

    extra_bits = bytes_needed * 8 - bits_needed
    if extra_bits > 0:
        num >>= extra_bits

    encoded = baseX_encode(num, alphabet)
    
    return encoded.rjust(out_length, alphabet[0])


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

