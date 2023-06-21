#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Merge HMF data into a set of Terra-style manifests
"""


import argparse
import json
import pandas as pd


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--metadata', help='HMF metadata.tsv',
                        required=True)
    parser.add_argument('--manifest', help='HMF manifest.json', 
                        required=True)
    parser.add_argument('-o', '--outdir', help='path to output directory',
                        default='./')
    args = parser.parse_args()

    # Load metadata (do not clean yet)
    md = pd.read_csv(args.metadata, sep='\t')
    import pdb; pdb.set_trace()
    
    # Parse .json and annotate WGS/RNAseq paths
    with open(args.manifest) as fin:
        for sdat in json.load(fin)['data']:
            import pdb; pdb.set_trace()


if __name__ == '__main__':
    main()

