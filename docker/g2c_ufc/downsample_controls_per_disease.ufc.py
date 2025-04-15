#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Adapted from Ryan L. Collins, Riaz Gillani, Jett Crowdis and the Van Allen Laboratory
# Copyright (c) 2023 Noah Fields and the Van Allen Laboratory
# Distributed under terms of the GPL-2.0 License (see LICENSE)
# Contact: Noah D. Fields <Noah_Fields@dfci.harvard.edu>

"""
Downsample control samples to ancestry-match individual pediatric cancer subtypes
"""


import argparse
import numpy as np
import pandas as pd
import random
import string
import networkx as nx
from scipy.stats import chisquare, fisher_exact


ascii_idx = {l : i + 1 for i, l in enumerate(string.ascii_lowercase)}

#### Noah's Version ####


#### End of Noah's Version ####

def melt_meta(meta, keepers):
    """
    Melt metadata pd.DataFrame into dict format expected by balance_ancestries()
    """

    def __build_entry(sdat):
        if sdat.disease == 'control' \
        and (sdat.proband == 'No' or sdat.isnull().proband):
            pheno = 'control'
        elif sdat.disease != 'control' \
        and (sdat.proband == 'Yes' or sdat.isnull().proband):
            pheno = sdat.disease
        else:
            pheno = None

        return {'pop' : sdat.ancestry_inferred_by_SVs,
                'pheno' : pheno,
                'case_control' : sdat.study_phase == 'case_control'}


    sample_dict = {}
    for idx, sdat in meta.loc[meta.sample_id.isin(keepers), :].iterrows():
         sample_dict[sdat.sample_id] = __build_entry(sdat)

    return sample_dict


def remove_sample(samples, sid):
    """
    Remove a sample from samples dict while accounting for relatives / "complete_trio" label
    """

    samples.pop(sid)

    for oid, values in samples.items():
        if 'trio_members' not in values.keys():
            continue
        if sid in values['trio_members']:
            samples[rel]['complete_trio'] = False

    return samples


def balance_ancestries(samples, require={}, seed=2023, abs_fold=True,
                       allow_drop_cases=True, case_weight=1,
                       return_keepers=False, verbose=False):
    """
    Downsample controls to balance ancestry distributions between cases and controls
    """

    # Subset samples based on requirements, if optioned
    samples_subset = samples.copy()
    for sid, vals in samples.items():
        for key, value in require.items():
            if isinstance(value, list):
                keep = vals[key] in value
            else:
                keep = vals[key] == value
            if not keep:
                samples_subset.pop(sid)
                break

    def __count_pops(samples):
        pops = {}
        for vals in samples.values():
            pop = vals['pop']
            if pop not in pops.keys():
                pops[pop] = {True : 0, False : 0}
            is_case = (vals['pheno'] != 'control')
            pops[pop][is_case] += 1
        return pd.DataFrame.from_dict(pops)

    def __compare_pops(counts):
        k_case = counts.loc[True, :].values
        k_ctrl = counts.loc[False, :].values
        ratio = sum(k_case) / sum(k_ctrl)
        x2 = chisquare(k_case, ratio * k_ctrl)
        return x2.pvalue

    def __get_fold(counts, abs_fold=True, biggest_step_always=True, case_weight=1):
        fold = counts.loc[False, :] / counts.loc[True, :]
        total_fold = counts.loc[False, :].sum() / counts.loc[True, :].sum()
        if abs_fold:
            if biggest_step_always:
                # When this option is enabled, we prioritize removing samples
                # from whichever population that will require the *fewest*
                # total samples removed to reach equilibrium with the other populations
                # Hence, the "biggest" step

                # Compute number of sample drops per population required 
                # to reach equilibrium with all other populations
                targets = pd.Series()
                for pop in fold.index:
                    other_counts = counts.loc[:, counts.columns != pop].sum(axis=1)
                    other_fold = other_counts.loc[False] / other_counts.loc[True]
                    n_case = counts.loc[True, pop]
                    n_ctrl = counts.loc[False, pop]
                    n_other_case = other_counts[True]
                    n_other_ctrl = other_counts[False]

                    if fold[pop] >= other_fold:
                        targets[pop] = n_ctrl - (n_case * other_fold)
                    else:
                        targets[pop] = case_weight * (n_case - (n_ctrl / other_fold))

                    # Test to ensure ratio is significantly different from other_fold
                    fisher_p = fisher_exact(np.array([other_counts, counts.loc[:, pop]]))[1]

                    # Set target to an absurdly large number if ratio is not 
                    # significantly different from all other pops (e.g., this pop)
                    # does not need to be pruned
                    if fisher_p >= 0.05:
                        targets[pop] = 10e6

                order = targets.sort_values(ascending=True).index

            else:
                order = np.log2(fold).abs().sort_values(ascending=False).index
        else:
            order = np.log2(fold).sort_values(ascending=False).index
        return fold[order], total_fold

    counts = __count_pops(samples_subset)
    if verbose:
        print('Samples prior to pruning:')
        print(counts)

    # Drop all samples from a population that has strictly zero cases or controls
    sids_to_prune = set()
    pop_is_quasi_empty = (counts == 0).any()
    for pop in pop_is_quasi_empty.index[pop_is_quasi_empty].tolist():
        for sid, vals in samples_subset.items():
            if vals['pop'] == pop:
                sids_to_prune.add(sid)
    for sid in sids_to_prune:
        sample_subset = remove_sample(samples_subset, sid)

    counts = counts.loc[:, ~pop_is_quasi_empty]
    pval = __compare_pops(counts)
    random.seed(seed)
    while pval < 0.01 or np.isnan(pval):
        fold, total_fold = __get_fold(counts, abs_fold=abs_fold, case_weight=case_weight)
        pop = fold.index[0]
        if fold[pop] >= total_fold or not allow_drop_cases:
            pop_sids = [sid for sid, vals in samples_subset.items() \
                        if vals['pop'] == pop and vals['pheno'] == 'control']
        else:
            pop_sids = [sid for sid, vals in samples_subset.items() \
                        if vals['pop'] == pop and vals['pheno'] != 'control']
        loser = random.choice(pop_sids)
        sids_to_prune.add(loser)
        sample_subset = remove_sample(samples_subset, loser)
        counts = __count_pops(samples_subset)
        pval = __compare_pops(counts)

    if verbose:
        print('Samples after pruning:')
        print(counts)

    if return_keepers:
        return samples_subset

    for sid in sids_to_prune:
        samples = remove_sample(samples, sid)

    return samples


def rebalance_cancer(samples, cancer, allow_drop_cases=False, case_weight=1):
    """
    Conduct ancestry matching separately for case:control and trio arms 
    """

    seed = int(''.join([str(ascii_idx[l]) for l in cancer]))

    if cancer == 'pancan':
        n_case_start = len([s for s, v in samples.items() if v['pheno'] in 'ewing neuroblastoma osteosarcoma'.split()])
    else:
        n_case_start = len([s for s, v in samples.items() if v['pheno'] == cancer])
    n_ctrl_start = len([s for s, v in samples.items() if v['pheno'] == 'control'])

    print('\n\nNow pruning ' + cancer + ':')
    case_keepers, ctrl_keepers = [], []

    def _pheno_match(svals, pheno):
        if cancer == 'pancan':
            return svals['pheno'] in 'ewing neuroblastoma osteosarcoma'.split()
        else:
            return svals['pheno'] == pheno

    if cancer == 'pancan':
        req_phenos = 'ewing neuroblastoma osteosarcoma control'.split()
    else:
        req_phenos = [cancer, 'control']        

    if len([v for v in samples.values() if _pheno_match(v, cancer) and v['case_control']]) > 0:

        # As of v2.5.2, use *all* controls for osteosarcoma cases
        # since there are no osteosarcoma trios
        if cancer == 'osteosarcoma':
            requirements = {'pheno' : req_phenos}
        else:
            requirements = {'pheno' : req_phenos, 'case_control' : True}

        print('\nDiscovery cohort summary:')
        case_control_samples = \
            balance_ancestries(samples.copy(), 
                               require=requirements, 
                               seed=seed, abs_fold=allow_drop_cases, 
                               allow_drop_cases=allow_drop_cases,
                               case_weight=case_weight, 
                               return_keepers=True, verbose=True)
        if cancer == 'pancan':
            case_keepers = [s for s, v in case_control_samples.items() if v ['pheno'] in 'ewing neuroblastoma osteosarcoma'.split()]
        else:
            case_keepers = [s for s, v in case_control_samples.items() if v ['pheno'] == cancer]
        ctrl_keepers += [s for s, v in case_control_samples.items() if v ['pheno'] == 'control']

    if len([v for v in samples.values() if _pheno_match(v, cancer) and not v['case_control']]) > 0:
        print('\nTrio/replication cohort summary:')
        trio_samples = \
            balance_ancestries(samples.copy(), 
                               require={'pheno' : req_phenos, 
                                        'case_control' : False}, 
                               seed=seed, abs_fold=allow_drop_cases, 
                               allow_drop_cases=allow_drop_cases, 
                               return_keepers=True, verbose=True)
        if cancer == 'pancan':
            case_keepers += [s for s, v in trio_samples.items() if v ['pheno'] in 'ewing neuroblastoma osteosarcoma'.split()]
        else:
            case_keepers += [s for s, v in trio_samples.items() if v ['pheno'] == cancer]
        ctrl_keepers += [s for s, v in trio_samples.items() if v ['pheno'] == 'control']
    
    n_case_end = len(case_keepers)
    n_ctrl_end = len(ctrl_keepers)

    print('\nRetained {:,} of {:,} cases for {}'.format(n_case_end, n_case_start, cancer))
    print('Retained {:,} of {:,} controls for {}'.format(n_ctrl_end, n_ctrl_start, cancer))

    return case_keepers, ctrl_keepers

def non_familial_set(samples, kinship_file):
    """
    Given a set of sample IDs and a kinship file (with columns #ID1 and ID2),
    return the subset of sample IDs that do NOT appear in either column.
    """
    kinship = pd.read_csv(kinship_file, delim_whitespace=True)

    # Collect all IDs that appear in the kinship pairs
    related_ids = set(kinship['#ID1']).union(set(kinship['ID2']))

    # Return the samples that are not related to anyone
    return samples - related_ids

def extract_family_units(kinship_file):
    """
    Given a kinship file with columns #ID1 and ID2,
    returns a list of sets, where each set represents a family unit
    (i.e., a connected component in the kinship graph).
    """
    kinship = pd.read_csv(kinship_file, delim_whitespace=True)
    
    G = nx.Graph()
    
    # Add edges from kinship file
    for _, row in kinship.iterrows():
        G.add_edge(row['#ID1'], row['ID2'])

    # Get connected components (family units)
    family_units = [set(component) for component in nx.connected_components(G)]
    
    return family_units

def dfs_independent_set(G, node_list, selected, best, cancer_status):
    # Base Case if the node list is empty
    if not node_list:
        selected_cases = {n for n in selected if cancer_status.get(n, 'control') != 'control'}
        best_cases = {n for n in best[0] if cancer_status.get(n, 'control') != 'control'}

        # Prefer more cancer cases; break ties with total size
        if (len(selected_cases) > len(best_cases)) or \
           (len(selected_cases) == len(best_cases) and len(selected) > len(best[0])):
            best[0] = selected.copy()
        return

    # Recursive Case where node list is not empty
    node = node_list[0]

    # Option 1: include node
    if all(neigh not in selected for neigh in G.neighbors(node)):
        selected.add(node)
        remaining = [n for n in node_list[1:] if n not in G.neighbors(node)]
        dfs_independent_set(G, remaining, selected, best, cancer_status)
        selected.remove(node)

    # Option 2: exclude node
    dfs_independent_set(G, node_list[1:], selected, best, cancer_status)



def maximal_non_related_subset_dfs(family, kinship_file, meta):
    """
    DFS-based search to find a maximal independent set within a family group.
    """
    # Load kinship
    kinship = pd.read_csv(kinship_file, delim_whitespace=True)
    edges = kinship[
        kinship['#ID1'].isin(family) & kinship['ID2'].isin(family)
    ][['#ID1', 'ID2']].values

    # Build graph
    G = nx.Graph()
    G.add_nodes_from(family)
    G.add_edges_from(edges)

    # Cancer info
    cancer_status = meta.set_index('original_id').loc[family]['cancer'].to_dict()
    nx.set_node_attributes(G, cancer_status, name='cancer')

    # Sort nodes to prioritize cases first
    node_list = sorted(
        G.nodes, key=lambda n: (cancer_status.get(n, 'control') == 'control', G.degree[n])
    )

    # DFS setup
    selected = set()
    best = [set()]  # use list as mutable container to hold best set
    dfs_independent_set(G, node_list, selected, best, cancer_status)

    return best[0]


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--metadata', required=True, help='sample metadata .tsv')
    parser.add_argument('--sample-list', required=True, help='list of samples to keep compared to larger G2C study')
    parser.add_argument('--cancer-subtype' required=True, help='cancer subtype being analyzed')
    parser.add_argument('--pca', required = True, help='.eigenvec file')
    parser.add_argument('--kinship', required = False, help='plink king file represtenting sparse matrix of kinship')
    parser.add_argument('--keep-samples', required = False, help='list of IDs to retain during ' +
                        'ancestry matching [default: keep all samples]')
    parser.add_argument('--outfile', required=True, help='output.tsv')
    args = parser.parse_args()

    # Load sample metadata
    meta = pd.read_csv(args.metadata, sep='\t',index=False)
    pca = pd.read_csv(args.pca,sep='\t',index=False)

    # Load list of samples to keep
    with open(args.sample_list) as f:
        samples = set(patient.strip() for patient in f)

    # Filter to just cases in our study as well as cases in the specific subtype
    meta = meta[meta['original_id'].isin(samples)]
    meta = meta[meta['cancer'] != "unknown"]
    if args.caner_subtype != "pancancer":
      meta = meta[meta['cancer'].str.contains(f"control|{args.cancer_subtype}")]

    ## Grab maximally unrelated set; enriching for cases ##
    # Grab all cases not involved in a family
    non_familial_set = non_familial_set(samples=set(meta['original_id']),kinship_file=args.kinship_file)
    family_units = extract_family_units(kinship_file=args.kinship_file)

    familial_set = set()
    for family in families:
        subset = maximal_non_related_subset_dfs(family, args.kinship_file, meta)
        familial_set.add(subset)
        print(subset)

    # Filter our data to our maximal unrelated set of individuals
    meta = meta[meta['original_id'].isin(non_familial_set.union(familial_set))]
    
    # Write to outfile
    meta.to_csv(args.outfile, sep='\t', index=False, na_rep='NA')

if __name__ == '__main__':
    main()