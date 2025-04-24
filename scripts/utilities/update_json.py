#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2024-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Updates one .json with values from a second .json, and formats to be more readable
"""


import argparse
import collections.abc
import json
from sys import stdin, stdout


def update(d, u):
    """
    Union-style dictionary update for nested dictionaries
    Taken from https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth
    """
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = update(d.get(k, {}), v)
        else:
            d[k] = v
    return d


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-json', help='input .json', default='stdin')
    parser.add_argument('-u', '--update-json', help='update .json', required=True)
    parser.add_argument('-o', '--output-json', help='output .json', default='stdout')
    args = parser.parse_args()

    # Open and read contents of input json
    if args.input_json in '- stdin /dev/stdin'.split():
        data = json.load(stdin)
    else:
        with open(args.input_json) as jin:
            data = json.load(jin)

    # Read contents of update json and update input json accordingly
    with open(args.update_json) as jup:
        updates = json.load(jup)
        data = update(data, updates)

    # Format updated json and print to --output-json
    if args.output_json in '- stdout /dev/stdout'.split():
        jout = stdout
    else:
        jout = open(args.output_json, 'w')
    json.dump(data, jout, indent=4, ensure_ascii=True, sort_keys=True)


if __name__ == '__main__':
    main()

