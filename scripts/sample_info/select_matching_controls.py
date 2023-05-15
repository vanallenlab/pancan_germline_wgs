#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Genomic Architecture of Human Cancers
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Randomly select covariate-matched controls from one cohort to match a target cohort of cases
"""


import argparse
import numpy as np
import pandas as pd
import random
from datetime import datetime
from itertools import product
from os import getcwd
from scipy.spatial.distance import pdist, squareform


def load_meta(meta_in, categorical_features, keep_columns=None):
    """
    Load & clean sample metadata
    """

    # Load data as pd.DataFrame & set sample ID as index
    meta = pd.read_csv(meta_in, sep="\t")
    meta.set_index(meta.columns.tolist()[0], drop=True, inplace=True)

    # Impute missing values for all non-categorical features as medians
    numeric_features = list(set(meta.columns.tolist()).difference(set(categorical_features)))
    numeric_fills = meta[numeric_features].median().to_dict()
    meta.fillna(numeric_fills, axis=0, inplace=True)

    # Impute missing values for categorical features as modes
    categorical_features = list(set(categorical_features).intersection(set(meta.columns.tolist())))
    categorical_fills = meta[categorical_features].mode().to_dict()
    meta.fillna(categorical_fills, axis=0, inplace=True)
    
    # Only return columns specified by keep_columns
    if keep_columns is not None:
        meta = meta.loc[:, keep_columns]
    return meta


def subset_samples(df, kv_pairs, drop_keys=True):
    """
    Subset a samples dataframe based on kv_pairs
    """

    for key, value in kv_pairs.items():
        df = df.loc[df[key] == value, :]

    if drop_keys:
        df.drop(kv_pairs.keys(), axis=1, inplace=True)

    return df


def match_controls(cases, controls, kv_pairs, numeric_features, 
                   ratio=1, seed=2023, quiet=True):
    """
    Select the best-matching N controls for each case after subsetting on kv_pairs
    """

    # Subset cases and controls based on kv_pairs before merging dataframes
    cases_sub = subset_samples(cases, kv_pairs)
    controls_sub = subset_samples(controls, kv_pairs)
    df = pd.concat([cases_sub, controls_sub], axis=0)

    # Report a cautionary warning if there are insufficient controls
    n_cases = cases_sub.shape[0]
    n_controls = controls_sub.shape[0]
    if not quiet and ratio * n_cases > n_controls:
        warn = '\nWarning: insufficient controls (N={:,}) to achieve {}:1 ratio ' + \
                'vs. cases (N={:,}) for the following condition:\n  {}'
        print(warn.format(n_controls, ratio, n_cases, kv_pairs))

    # Compute distance matrix between all cases and controls
    # If any numeric features are present, standard normalize all variables before
    # computing Euclidean distance of feature vectors
    if len(df.shape) == 2 and df.shape[-1] > 0:

        df = (df - df.mean()) / df.std()

        # Compute distance matrix between all samples
        dist = pd.DataFrame(squareform(pdist(df)), 
                            columns=df.index, index=df.index)
        dist = dist.loc[cases_sub.index, controls_sub.index]

    # If no numeric features are present, arbitrarily set all distances to the same value (0)
    else:
        dist = pd.DataFrame(0, index=cases_sub.index, columns=controls_sub.index)

    # Randomly shuffle row and column orders to ensure random tie breaks
    random.seed(seed)
    dist = dist.loc[random.sample(dist.index.tolist(), dist.shape[0]), 
                    random.sample(dist.columns.tolist(), dist.shape[1])]

    # Matching algorithm
    keepers = {case_id : set() for case_id in cases_sub.index}
    # Iterate until all cases have exactly N (ratio) controls
    # OR no controls remain
    while all([k > 0 for k in dist.shape]):

        # Find case with worst best match (i.e., largest row-wise min)
        # If multiple rows are tied, arbitrarily choose the first one
        best_dist = dist.min(axis=1)
        case_id = best_dist.index[best_dist == np.nanmax(best_dist)][0]
        
        # Find the best control(s) for the worst case
        # If multiple columns are tied, arbitrarily choose the first one
        control_dists = dist.loc[case_id, :]
        control_id = control_dists.index[control_dists == np.nanmin(control_dists)][0]

        # Add new case-control pair to keepers and drop control from dist matrix
        # Also drop case if its quota has been reached
        keepers[case_id].add(control_id)
        dist.drop(control_id, axis=1, inplace=True)
        if len(keepers[case_id]) == ratio:
            dist.drop(case_id, axis=0, inplace=True)

    return keepers


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('case_meta', help='input case metadata .tsv')
    parser.add_argument('control_meta', help='input control metadata .tsv')
    parser.add_argument('outfile', help='output list of subsetted sample IDs')
    parser.add_argument('--annotate-outfile', action='store_true',
                        help='Add second column to outfile specifying the case ' +
                        'sample ID to which each control was matched.')
    parser.add_argument('-n', '--n-per-case', default=1, type=int,
                        help='Number of controls to match per case.')
    parser.add_argument('--categorical-variable', action='extend', nargs='?',
                        help='Specify custom variable to be categorically matched. ' +
                        'Can be provided multiple times.')
    parser.add_argument('--seed', type=int, default=2023, help='Random seed')
    parser.add_argument('--quiet', action='store_true', help='Suppress logs.')
    parser.add_argument
    args = parser.parse_args()

    # Start log unless disabled
    if not args.quiet:
        log_msg = '\n*** Log for select_matching_controls.py ***\n\n'
        log_msg += 'Execution directory: {}\n'.format(getcwd())
        log_msg += 'Execution date and time: {}\n'.format(datetime.now().strftime("%m/%d/%Y at %H:%M:%S"))
        log_msg += 'Arguments:'
        print(log_msg)
        for key, value in vars(args).items():
            print('  - {} : {}'.format(key, value))
        print('')

    # Build list of all categorical variables
    categorical_features = 'population ancestry race sex gender'.split()
    if args.categorical_variable is not None:
        categorical_features += args.categorical_variable
    categorical_features = sorted([f.lower() for f in set(categorical_features)])

    # Load all case metadata
    cases = load_meta(args.case_meta, categorical_features)
    if not args.quiet:
        log_msg = 'Loaded {:,} features for {:,} cases from {}'
        print(log_msg.format(*cases.shape[::-1], args.case_meta))

    # Load all control metadata
    controls = load_meta(args.control_meta, categorical_features, cases.columns)
    if not args.quiet:
        log_msg = 'Loaded {:,} features for {:,} controls from {}'
        print(log_msg.format(*controls.shape[::-1], args.case_meta))
        print('')

    # Finalize features to use during matching
    categorical_features = list(set(cases.columns.tolist()).intersection(set(categorical_features)))
    numeric_features = list(set(cases.columns.tolist()).difference(set(categorical_features)))
    if not args.quiet:
        print('Matching on the following features:')
        if len(categorical_features) > 0:
            print('  - Categorical: {}'.format(', '.join(categorical_features)))
        if len(numeric_features) > 0:
            print('  - Numeric: {}'.format(', '.join(numeric_features)))

    # Build list of observed values for each categorical feature
    categorical_values = {f : cases[f].unique().tolist() for f in categorical_features}

    # Perform matching
    keepers = {}
    for vals in product(*categorical_values.values()):
        kv_pairs = {k : v for k, v in zip(categorical_values.keys(), vals)}
        new_keepers = match_controls(cases, controls, kv_pairs, numeric_features, 
                                     args.n_per_case, args.seed, args.quiet)
        keepers.update(new_keepers)

    # Write out keepers to file
    with open(args.outfile, 'w') as fout:
        for case, controls in keepers.items():
            for control in controls:
                if args.annotate_outfile:
                    fout.write('{}\t{}\n'.format(control, case))
                else:
                    fout.write(control + '\n')


if __name__ == '__main__':
    main()

