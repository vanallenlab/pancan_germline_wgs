#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Genomic Architecture of Human Cancers
# Copyright (c) 2023-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Perform sample batching and outlier exclusion for a single cohort before GATK-SV pipeline
"""


import argparse
import json
import numpy as np
import operator as op
import pandas as pd
from datetime import datetime
from itertools import product
from os.path import basename
from scipy.stats import median_abs_deviation as mad
from sys import stdout


# Set global variables used in a variety of scopes
qc_log_fmt = ' - {:,} ({:.1f}%) samples failed due to {}\n'
eqs = {'<' : op.lt, '<=' : op.le, '==' : op.eq, 
       '!=' : op.ne, '>=' : op.ge, '>' : op.gt}
maths = {'q1' : lambda x: np.nanquantile(x, q=0.25), 'mean' : np.nanmean, 
         'median' : np.nanmedian, 'q3' : lambda x: np.nanquantile(x, q=0.75), 
         'sd' : np.nanstd, 'mad' : mad}


def print_log_header(args, log):
    """
    Print header for logfile
    """

    title = '### GATK-SV Sample Batching Tool ###'
    spacer = '#' * len(title)
    log.write('\n'.join([spacer, title, spacer]) + '\n')
    log.write('\n' + datetime.now().strftime('%B %d, %Y at %H:%M\n'))
    log.write('\nRuntime arguments:\n')
    for arg in vars(args):
        log.write(' - {} : {}\n'.format(arg, getattr(args, arg)))


def load_metadata(meta_in, log, quiet=False):
    """
    Read sample metadata
    """

    md = pd.read_csv(meta_in, sep='\t', low_memory=False)
    md.set_index(md.columns[0], drop=True, inplace=True)
    
    if not quiet:
        msg = '\nLoaded {:,} samples and {:,} features from {}\n'
        log.write(msg.format(*md.shape, basename(meta_in)) + '\n')
    
    return md


def parse_cutoff(cut_str, values):
    """
    Parse a string expression for a filter cutoff
    """

    # Sequentially replace all special strings with their numeric equivalents
    for key, fx in maths.items():
        cut_str = cut_str.replace(key, str(fx(values)))

    # Attempt to evaluate string as mathematical expression
    return eval(cut_str)


def run_qc(md, jsonfile, log, quiet=False):
    """
    Apply filters specified in json_in to pd.Dataframe md (sample metadata)
    Returns a set of samples failing any filter
    If quiet = False, will also report number of samples failing each filter
    """

    fails = set()
    for filt_name, filt_info in json.load(jsonfile).items():
        feature = filt_info.get("feature")
        # Check if feature is present in md columns
        if feature not in md.columns:
            msg = 'Feature "{}" is not present in sample metadata header ' + \
                  'but was specified by {} in {}. Exiting.'
            exit(msg.format(feature, filt_name, jsonfile.name))

        # Parse equality
        eq_fx = eqs.get(filt_info['equality'])
        if eq_fx is None:
            msg = 'Equality "{}" was specified for feature "{}" by {} in {} but could ' + \
                  'not be interpreted. Exiting.'
            exit(msg.format(filt_info['equality'], feature, filt_name, jsonfile.name))

        # Parse cutoff
        cutoff = parse_cutoff(filt_info['value'].lower(), md[feature])

        # Apply filter
        new_fails = set(md.index[md[feature].apply(lambda v: eq_fx(v, cutoff))].tolist())
        fails.update(new_fails)
        if not quiet:
            fail_reason = '{} ({} {} {:.2f} [{}])'.format(filt_name, filt_info['feature'],
                                                          filt_info['equality'], cutoff, 
                                                          filt_info['value'])
            log.write(qc_log_fmt.format(len(new_fails), 100 * len(new_fails) / md.shape[0],
                                        fail_reason))

    return fails


def filter_samples(md, criteria):
    """
    Extract the IDs for samples meeting all criteria specified as key : value pairs
    """

    # Always restrict to QC-passing samples
    criteria = dict(criteria)
    for key in 'global_qc_pass batch_qc_pass'.split():
        criteria[key] = True

    keepers = pd.Series([True] * md.shape[0], index=md.index)
    for key, value in criteria.items():
        if key in md.columns:
            keepers = keepers & (md[key] == value)
    
    return list(keepers.index[keepers])


def summarize_batches(array, log):
    """
    Summarize distribution of batch sizes
    """

    k = array.value_counts().value_counts()
    for size, n in k.sort_index().to_dict().items():
        log.write(' - {:,} batches @ {:,} samples\n'.format(n, size))
    log.write('\n')


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('sample_metadata', help='Sample metadata .tsv')
    parser.add_argument('--match-on', action='append', type=str, default=[],
                        help='Column name from sample_metadata to be treated ' +
                        'as a categorical stratum when batching such that all ' +
                        'resulting batches have ~the same number of samples ' +
                        'from each strata. Will be applied before --batch-by.')
    parser.add_argument('--batch-by', action='append', type=str, default=[],
                        help='Column name from sample_metadata to use when ' +
                        'dividing samples into batches. Note: order matters! ' +
                        'Samples will be batched based on the order of features ' +
                        'provided to --batch-by.')
    parser.add_argument('--global-qc-cutoffs', help='Optional .json specifying ' +
                        'QC cutoffs for global sample exclusion prior to batching.')
    parser.add_argument('--custom-qc-fail-samples', help='Optional list of samples ' +
                        'to be treated as failing global QC and excluded from ' +
                        'batching.')
    parser.add_argument('--batch-qc-cutoffs', help='Optional .json specifying ' +
                        'QC cutoffs for batch-specific sample exclusion.')
    parser.add_argument('--batch-size', type=int, default=500,
                        help='ideal batch size')
    parser.add_argument('--prefix', default='wgs_batch',
                        help='prefix for all batches')
    parser.add_argument('-o', '--outfile', required=True, help='path to revised ' +
                        'metadata .tsv output with QC and batch labels')
    parser.add_argument('--batch-names-tsv', help='path to Terra-style .tsv file ' + 
                        'output with names of all batches')
    parser.add_argument('--batch-membership-tsv', help='prefix for Terra-style .tsv ' +
                        'file with IDs of samples per batch')
    parser.add_argument('-l', '--logfile', default='stdout', help='Path to ' +
                        'diagnostic log file. Suppressed by --quiet.')
    parser.add_argument('-q', '--quiet', action='store_true', help='Suppress ' +
                        'verbose logging.')
    args = parser.parse_args()

    # Check that at least one batching variable was specified
    if len(args.batch_by) == 0:
        exit('At least one value of --batch-by must be specified. Exiting.')

    # Sanity check all variables are present in sample_metadata
    with open(args.sample_metadata) as fin:
        colnames = [h.replace('#', '') for h in fin.readline().rstrip().split('\t')]
    matchon_missing = set(args.match_on).difference(set(colnames))
    batchby_missing = set(args.batch_by).difference(set(colnames))
    missing_colname_msg = 'The following features were passed to {} but could ' + \
                           'be located in the input sample metadata .tsv: {}'
    if len(matchon_missing) > 0:
        err(missing_colname_msg.format('--match-on', ', '.join(matchon_missing)))
    if len(batchby_missing) > 0:
        err(missing_colname_msg.format('--batch-by', ', '.join(batchby_missing)))

    # Print log title
    if args.logfile in 'stdout /dev/stdout -'.split():
        log = stdout
    else:
        log = open(args.logfile, 'w')
    if not args.quiet:
        print_log_header(args, log)

    # Load sample metadata
    md = load_metadata(args.sample_metadata, log, args.quiet)
    qc_args = 'global_qc_cutoffs custom_qc_fail_samples batch_qc_cutoffs'.split()
    if any(getattr(args, x) is None for x in qc_args):
        md['preqc_batch_assignment'] = [None] * md.shape[0]

    # Global sample QC, if optioned
    if not all(getattr(args, x) is None for x in qc_args[:2]):
        md['global_qc_pass'] = [True] * md.shape[0]
        global_fails = set()
        if not args.quiet:
            log.write('Performing global QC before batching. Summary:\n')
        if args.global_qc_cutoffs is not None:
            with open(args.global_qc_cutoffs) as j_in:
                global_qc_fails = run_qc(md, j_in, log, args.quiet)
                global_fails.update(global_qc_fails)
        if args.custom_qc_fail_samples is not None:
            with open(args.custom_qc_fail_samples) as fin:
                hard_fails = set([s.rstrip() for s in fin.readlines()])
                global_fails.update(hard_fails)
                n_hard_fails = len(hard_fails)
                if not args.quiet:
                    qc_log_vals = [n_hard_fails, 100 * n_hard_fails / md.shape[0],
                                   'being specified in --custom-qc-fail-samples']
                    log.write(qc_log_fmt.format(*qc_log_vals))
        if not args.quiet:
            qc_log_vals = [len(global_fails), 100 * len(global_fails) / md.shape[0],
                           'any of the above reasons']
            log.write(qc_log_fmt.format(*qc_log_vals))
            msg = 'Retained a total of {:,} samples after global QC.\n\n'
            log.write(msg.format(md.shape[0] - len(global_fails)))
        md.loc[list(global_fails), 'global_qc_pass'] = False

    else:
        if not args.quiet:
            log.write('No --global-qc-cutoffs supplied; skipping pre-batching sample QC\n\n')

    # Enumerate matching strata
    strata_values = [list(md[f].unique()) for f in args.match_on]
    strata = {k : set() for k in tuple(product(*strata_values))}
    for ks in strata.keys():
        strata[ks] = filter_samples(md, zip(args.match_on, ks))
    if not args.quiet and len(strata_values) > 0:
        msg = 'Identified {:,} strata based on unique combinations of --match-on:\n'
        log.write(msg.format(len(strata)))
        msg = ' - {} : {:,} samples\n'
        for ks, samps in strata.items():
            log.write(msg.format(' '.join(ks), len(samps)))
        log.write('Batching each stratum in parallel before pooling across strata.\n\n')

    # Estimate batching partitions
    bb_numeric = set([f for f in args.batch_by if pd.api.types.is_numeric_dtype(md[f])])
    bb_categ = {f : list(md[f].unique()) for f in set(args.batch_by).difference(bb_numeric)}
    # Since the number of splits per numeric feature is adjustable but the number
    # of splits per categorical feature is not, we first need to count the number of
    # samples for each unique combination of categorical features before deciding on
    # number of splits per numeric feature
    batch_trees = {}
    for categ_values in product(*bb_categ.values()):
        branches = {k : [0] for k in args.batch_by}
        for key, value in zip(bb_categ, categ_values):
            branches[key] = [value]
        n_samples = len(filter_samples(md, zip(bb_categ.keys(), categ_values)))
        n_splits = np.prod([len(v) for v in branches.values()])
        prev_batch_size = current_batch_size = n_samples / n_splits
        while current_batch_size > args.batch_size:
            # Using a round-robin strategy, find the highest-priority numeric
            # --split-by feature that has the fewest number of splits, and add
            # one split to that feature
            current_min_splits = np.nanmin([len(branches[k]) for k in bb_numeric])
            for k in args.batch_by:
                if k in bb_numeric and len(branches.get(k, [])) == current_min_splits:
                    branches[k].append(max(branches[k]) + 1)
                    prev_batch_size = current_batch_size
                    n_splits = np.prod([len(v) for v in branches.values()])
                    current_batch_size = n_samples / n_splits
                    break
        # Once optimum number of splits has been straddled, select either the current
        # or previous step based on whichever results in a batch size closer to
        # the specified target --batch-size
        d_current = abs(args.batch_size - current_batch_size)
        d_prev = abs(args.batch_size - prev_batch_size)
        if d_current > d_prev:
            branches[k] = branches[k][:-1]
        batch_trees[categ_values] = branches
    
    # Perform batching on parallel strata from above
    batch_name_template = args.prefix
    for bb in args.batch_by:
        if bb in bb_numeric:
            batch_name_template += '|' + bb + '_{}'
        else:
            batch_name_template += '|{}'
    for strat, strat_ids in strata.items():
        for cv, branches in batch_trees.items():
            samples = filter_samples(md, zip(bb_categ.keys(), cv))
            samples = set(samples).intersection(strat_ids)
            for combo in product(*branches.values()):
                bnvals = []
                for v in combo:
                    if str(v).isnumeric():
                        bnvals.append(v + 1)
                    else:
                        bnvals.append(v)
                batch_name = batch_name_template.format(*bnvals)
                md_sub = md.loc[list(samples), args.batch_by].copy()
                for key, value in zip(branches.keys(), combo):
                    if key in bb_numeric:
                        max_quant = len(branches[key])
                        md_sub[key] = np.floor(max_quant * md_sub[key].rank(method='first') / (md_sub.shape[0] + 1))
                    md_sub = md_sub.loc[md_sub[key] == value, :]
                md.loc[md_sub.index.tolist(), 'preqc_batch_assignment'] = batch_name

    # Report batching summary prior to QC
    if not args.quiet:
        log.write('Batching finished. Batch membership summary:\n')
        summarize_batches(md['preqc_batch_assignment'], log)

    # Batch-specific sample QC, if optioned
    md['final_batch_assignment'] = md['preqc_batch_assignment'].copy()
    if args.batch_qc_cutoffs is None:
        if not args.quiet:
            log.write('No --batch-qc-cutoffs supplied; skipping batch-specific sample QC\n\n')
    else:
        md['batch_qc_pass'] = md['global_qc_pass'].copy()
        batch_qc_fails = set()
        batch_names = md.loc[~md['preqc_batch_assignment'].isna(), 
                             'preqc_batch_assignment'].unique()
        if not args.quiet:
            log.write('Performing batch-specific QC. Summary per batch:\n\n')
        for batch in sorted(batch_names):
            batch_md = md.loc[md['preqc_batch_assignment'] == batch, :]
            if not args.quiet:
                log.write('QC summary for batch "{}":\n'.format(batch))
            with open(args.batch_qc_cutoffs) as j_in:
                new_batch_qc_fails = run_qc(batch_md, j_in, log, args.quiet)
            if not args.quiet:
                qc_log_vals = [len(new_batch_qc_fails), 
                               100 * len(new_batch_qc_fails) / batch_md.shape[0],
                               'any of the above reasons']
                log.write(qc_log_fmt.format(*qc_log_vals))
                msg = 'Retained a total of {:,} samples in batch "{}" after QC.\n\n'
                log.write(msg.format(batch_md.shape[0] - len(new_batch_qc_fails), batch))
            batch_qc_fails.update(new_batch_qc_fails)
        md.loc[list(batch_qc_fails), 'batch_qc_pass'] = False
        md.loc[list(batch_qc_fails), 'final_batch_assignment'] = None
        if not args.quiet:
            n_kept = md['batch_qc_pass'].sum()
            msg = 'Batch-specific QC complete. Retained a total of {:,} ' + \
                  'samples across all batches.\n\n'
            log.write(msg.format(n_kept))
            log.write('Final batch membership membership summary after QC:\n')
            summarize_batches(md['final_batch_assignment'], log)

    # Write updated metadata to --outfile
    md.insert(0, colnames[0], md.index)
    md.to_csv(args.outfile, sep='\t', index=False, header=True, na_rep='NA')

    # Write Terra-formatted .tsvs, if optioned
    if args.batch_names_tsv is not None:
        with open(args.batch_names_tsv, 'w') as fout:
            fout.write('entity:sample_set_id\n')
            for batch in sorted(md.loc[md.final_batch_assignment.notna(), 
                                       'final_batch_assignment'].unique()):
                fout.write(batch + '\n')
    if args.batch_membership_tsv is not None:
        with open(args.batch_membership_tsv, 'w') as fout:
            fout.write('membership:sample_set_id\tsample\n')
            for batch in sorted(md.loc[md.final_batch_assignment.notna(), 
                                       'final_batch_assignment'].unique()):
                for sid in md.index[md['final_batch_assignment'] == batch].tolist():
                    fout.write('{}\t{}\n'.format(batch, sid))

    if not args.quiet:
        log.write('Finished.\n')
    log.close()


if __name__ == '__main__':
    main()
