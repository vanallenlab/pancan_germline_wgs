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


#ascii_idx = {l : i + 1 for i, l in enumerate(string.ascii_lowercase)}

#### Noah's Version ####

def extract_non_familial_set(samples, kinship_file):
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

def initial_filter(df,min_age=0,max_age=200,sex_karyotypes={"XX","XY"},apparent_aneuploidies={0}):

    cohort_set = {'aou','ceph','mesa','ufc'}
    df = df[df['cohort'].isin(cohort_set)]
    df = df[df['sex_karyotype'].isin(sex_karyotypes)]
    df = df[df['apparent_aneuploidies'].isin(apparent_aneuploidies)]
    df = df[df['age'] >= min_age]
    df = df[df['age'] <= max_age]
    df.to_csv("blah.tsv",sep='\t',index=False)
    return df

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
    Only includes individuals present in the metadata.
    """

    # Filter family members to only those present in metadata
    available_ids = set(meta['original_id'])
    family = family & available_ids  # intersection

    # Load kinship
    kinship = pd.read_csv(kinship_file, delim_whitespace=True)
    edges = kinship[
        kinship['#ID1'].isin(family) & kinship['ID2'].isin(family)
    ][['#ID1', 'ID2']].values


    if not family:
        return set()  # nothing to do if no valid individuals

    # Build graph
    G = nx.Graph()
    G.add_nodes_from(family)
    G.add_edges_from(edges)

    # Cancer info
    cancer_status = meta.set_index('original_id').loc[list(family)]['cancer'].to_dict()
    #cancer_status = meta.set_index('original_id').loc[family]['cancer'].to_dict()
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
    parser.add_argument('--cancer-subtype', required=True, help='cancer subtype being analyzed')
    parser.add_argument('--pca', required = True, help='.eigenvec file')
    parser.add_argument('--kinship', required = False, help='plink king file represtenting sparse matrix of kinship')
    parser.add_argument('--keep-samples', required = False, help='list of IDs to retain during ' +
                        'ancestry matching [default: keep all samples]')
    parser.add_argument('--min-age',required = False, help='minimum age for the study')
    parser.add_argument('--max-age',required = False, help='maximum age for the study')
    parser.add_argument('--apparent_aneuploidies',required = False, help='allowed ploidies')
    parser.add_argument('--sex-karyotypes',required = False, help='allowed sex ploidies')
    parser.add_argument('--outfile', required=True, help='output.tsv')
    parser.add_argument('--exclude-samples',required=False, help='Samples to Exclude')
    parser.add_argument('--log-file', required=False, help='Path to log file')
    args = parser.parse_args()

    # Make the log file
    if not args.log_file:
        args.log_file = f"{args.cancer_subtype}.cohort.log"

    with open(args.log_file, "w") as f:
        f.write("Size\tNum_Filtered\tPercent_Filtered\tExclusion_Criteria\n")

    # Load sample metadata
    meta = pd.read_csv(args.metadata, sep='\t',index_col=False)
    pca = pd.read_csv(args.pca,sep='\t',index_col=False)
    
    # Load list of samples to keep
    with open(args.sample_list) as f2:
        samples = set(patient.strip() for patient in f2)

    # Load list of samples to exclude
    if args.exclude_samples:
        with open(args.exclude_samples) as f3:
            exclude_samples = set(patient.strip() for patient in f3)
        samples = samples - exclude_samples

    # Filter to just cases in our study as well as cases in the specific subtype
    meta = meta[meta['original_id'].astype(str).str.strip().isin(samples)]
    sample_size1 = len(meta)

    with open(args.log_file, "a") as f:
        f.write(f"{sample_size1}\t0\t0\tInitial Samples\n")

    # Filter to only samples with known cancer status
    meta = meta[meta['cancer'] != "unknown"]
    sample_size2 = len(meta)
    
    with open(args.log_file, "a") as f:
        f.write(f"{sample_size2}\t{(sample_size1 - sample_size2)}\t{((sample_size1 - sample_size2)/sample_size1)}\tRemove samples with unknown cancer status.\n")


    if args.cancer_subtype != "pancancer":
      meta =meta[meta['cancer'].str.contains(f"control|{args.cancer_subtype}")]

    # Remove samples with irrelevant cancer diagnosis for this study
    sample_size3 = len(meta)

    with open(args.log_file, "a") as f:
        f.write(f"{sample_size3}\t{(sample_size2 - sample_size3)}\t{((sample_size2 - sample_size3)/sample_size2)}\tRemove samples that are not controls nor {args.cancer_subtype.replace('|',',')} diagnosis.\n") 

    # Do initial filtering of dataset
    if args.sex_karyotypes:
        study_sex_karyotypes = set(args.sex_karyotypes.split(','))
        meta = initial_filter(meta,sex_karyotypes=study_sex_karyotypes)
    else:
        meta = initial_filter(meta)

    # Remove samples with irrelevant cancer diagnosis for this study
    sample_size4 = len(meta)

    with open(args.log_file, "a") as f:
        f.write(f"{sample_size4}\t{(sample_size3 - sample_size4)}\t{((sample_size3 - sample_size4)/sample_size3)}\tRemove samples that are not {args.sex_karyotypes}.\n")

    ## Grab maximally unrelated set; enriching for cases ##
    # Grab all cases not involved in a family
    non_familial_set = extract_non_familial_set(samples=set(meta['original_id']),kinship_file=args.kinship)
    family_units = extract_family_units(kinship_file=args.kinship)

    familial_set = set()
    for family in family_units:
        subset = maximal_non_related_subset_dfs(family, args.kinship, meta)
        familial_set.update(subset)

    # Filter our data to our maximal unrelated set of individuals
    meta = meta[meta['original_id'].isin(non_familial_set.union(familial_set))]

    sample_size5 = len(meta)
    with open(args.log_file, "a") as f:
        f.write(f"{sample_size5}\t{(sample_size4 - sample_size5)}\t{((sample_size4 - sample_size5)/sample_size4)}\tExcluded {len(familial_set)} due to relatedness with other individuals.")

    
    ## Print Summary Statistics
    with open(args.log_file, 'a') as f:
        f.write("===== Summary Report =====\n\n")

        # Total number of individuals
        total_count = len(meta)
        f.write(f"Total individuals: {total_count}\n\n")

        # Count per cohort
        f.write("Counts by cohort:\n")
        f.write(meta['cohort'].value_counts().to_string())
        f.write("\n\n")

        # Count per cancer
        f.write("Counts by cancer:\n")
        f.write(meta['cancer'].value_counts().to_string())
        f.write("\n\n")

        # Count per sex_karyotype
        f.write("Counts by sex_karyotype:\n")
        f.write(meta['sex_karyotype'].value_counts().to_string())
        f.write("\n\n")

        # Mean age and BMI stratified by cancer type
        f.write("Mean age by cancer type:\n")
        means = meta.groupby('cancer')[['age']].mean().round(2)
        f.write(means.to_string())
        f.write("\n\n")

        # intake_qc_pop count per cancer
        f.write("Counts of intake_qc_pop per cancer:\n")
        intake_counts = meta.groupby('cancer')['intake_qc_pop'].value_counts().unstack(fill_value=0)
        f.write(intake_counts.to_string())
        f.write("\n\n")

        f.write("==========================\n\n")


    meta = meta.merge(pca, left_on='original_id', right_on='#IID', how='left')

    # Write to outfile
    meta.to_csv(args.outfile, sep='\t', index=False, na_rep='NA')

if __name__ == '__main__':
    main()
