#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Adapted from Ryan L. Collins, Riaz Gillani, Jett Crowdis and the Van Allen Laboratory
# Copyright (c) 2025 Noah Fields and the Van Allen Laboratory
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
import math
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

def initial_filter(df,min_age=0,max_age=200,sex_karyotypes={"XX","XY"},apparent_aneuploidies={0},cohorts="aou"):

    cohort_set = set(cohorts.split(','))
    df = df[df['cohort'].isin(cohort_set)]
    df = df[df['sex_karyotype'].isin(sex_karyotypes)]
    df = df[df['apparent_aneuploidies'].isin(apparent_aneuploidies)]
    df = df[df['age'] >= min_age]
    df = df[df['age'] <= max_age]
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

def match_controls_by_ancestry(meta: pd.DataFrame, max_controls_per_case: int = 3, seed: int = 42) -> pd.DataFrame:
    """
    Match controls to cases based on ancestry (intake_qc_pop), sex, and cohort if available.
    Try to balance age within each matched group.
    
    Parameters:
        meta: DataFrame with at least ['cancer', 'intake_qc_pop', 'sex', 'age', 'cohort']
        max_controls_per_case: Max number of controls to retain per case (default: 3)
        seed: Seed for reproducible sampling

    Returns:
        DataFrame of all cases and matched controls.
    """
    np.random.seed(seed)

    already_matched_control_ids = set()

    # Define control and case status
    is_control = meta['cancer'].str.lower().str.contains('control')
    is_case = ~meta['cancer'].isin(['control', 'unknown'])

    cases = meta[is_case].copy()
    controls = meta[is_control].copy()
    matched_controls = []

    # Iterate over ancestry + sex combinations in cases
    for (pop, sex), case_group in cases.groupby(['intake_qc_pop', 'sex_karyotype']):
        case_count = len(case_group)

        # To account for ancestry groups with disproportionately few cases, 
        # we cap the control-to-case ratio at 3:1 when such imbalances are present. 
        if (max_controls_per_case < 3):
            max_controls = case_count * 3
        else:
            max_controls = case_count * max_controls_per_case

        median_age = case_group['age'].median()

        # Find matching controls for ancestry and sex
        potential_controls = controls[
            (controls['intake_qc_pop'] == pop) &
            (controls['sex_karyotype'] == sex) &
            (~controls['original_id'].isin(already_matched_control_ids))
        ].copy()

        if potential_controls.empty:
            continue

        # Prioritize controls by age difference
        potential_controls['age_diff'] = (potential_controls['age'] - median_age).abs()
        potential_controls = potential_controls.drop_duplicates(subset='original_id')
        selected = potential_controls.sort_values('age_diff').head(max_controls)

        already_matched_control_ids.update(selected['original_id'])
        matched_controls.append(selected)

    # Combine everything
    matched_controls_df = pd.concat(matched_controls, ignore_index=True)
    final_meta = pd.concat([cases, matched_controls_df], ignore_index=True)

    return final_meta


def get_min_case_control_ratio(meta: pd.DataFrame):
    """
    Returns the case:control ratio for each ancestry group and identifies the smallest ratio.
    
    Parameters:
        meta: DataFrame with 'cancer' and 'intake_qc_pop' columns
    
    Returns:
        A Series with case/control ratios per ancestry, and a printed note about the smallest one.
    """
    # Define cases and controls
    meta = meta.copy()
    meta['is_case'] = ~meta['cancer'].str.lower().isin(['control', 'unknown'])

    # Group by ancestry
    grouped = meta.groupby('intake_qc_pop')['is_case']

    # Count cases and controls
    counts = grouped.agg(['sum', 'count'])  # sum = cases, count - sum = controls
    counts['controls'] = counts['count'] - counts['sum']
    counts['control_case_ratio'] =  counts['controls'] / counts['sum']

    # Handle div-by-zero (if any population has 0 controls)
    counts['control_case_ratio'] = counts['control_case_ratio'].replace([float('inf')], float('nan'))

    # Drop intermediate columns for clarity
    ratio_series = counts['control_case_ratio']
    
    # Print info on the ancestry with the lowest ratio
    min_pop = ratio_series.idxmin()
    min_val = ratio_series.min()
    print(f"Smallest control:case ratio is {min_val:.3f} in population: {min_pop}")
    print(f"Ratio Series: {ratio_series}")
    return min_val

def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--metadata', required=True, help='sample metadata .tsv')
    parser.add_argument('--phenotype-data', required=True, help='data describing phenotypes')
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
    parser.add_argument('--cohorts',required=True, help='Comma delimited string denoting which cohorts to include')
    parser.add_argument('--use-original-dx', required=True, default=False, help='Option to use original_dx to get more niche subtypes.')
    args = parser.parse_args()

    # Make the log file
    if not args.log_file:
        args.log_file = f"{args.cancer_subtype}.cohort.log"

    with open(args.log_file, "w") as f:
        f.write("Size\tNum_Filtered\tPercent_Filtered\tExclusion_Criteria\n")
        

    # Load sample metadata
    meta = pd.read_csv(args.metadata, sep='\t',index_col=False)
    phenotype_data = pd.read_csv(args.phenotype_data,sep='\t',index_col=False)
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

    # Merge the metadata with phenotype data
    phenotype_data = phenotype_data.drop(columns=['cancer','age'], errors='ignore')
    meta = meta.merge(phenotype_data, left_on = "original_id", right_on = "Sample")
    sample_size1 = len(meta)

    with open(args.log_file, "a") as f:
        f.write(f"{sample_size1}\t0\t0\tInitial Samples\n")

    # Filter to only samples with known cancer status
    meta = meta[meta['cancer'] != "unknown"]
    sample_size2 = len(meta)
    
    with open(args.log_file, "a") as f:
        f.write(f"{sample_size2}\t{(sample_size1 - sample_size2)}\t{round(((sample_size1 - sample_size2)/sample_size1),3) * 100}\tRemove samples with unknown cancer status.\n")


    # Filter to Cancer Subtypes of Interest
    print("Args Variables")
    print(vars(args))
    if args.cancer_subtype != "pancancer":
        if not args.use_original_dx:
            meta =meta[meta['cancer'].str.contains(f"control|{args.cancer_subtype}")]
        else:
            meta =meta[(meta['cancer'] == "control") | (meta['original_dx'].str.contains(args.cancer_subtype))]

    # Remove samples with irrelevant cancer diagnosis for this study
    sample_size3 = len(meta)

    with open(args.log_file, "a") as f:
        f.write(f"{sample_size3}\t{(sample_size2 - sample_size3)}\t{round(((sample_size2 - sample_size3)/sample_size2),3) * 100}\tRemove samples that are not controls nor {args.cancer_subtype.replace('|',',')} diagnosis.\n") 

    # Do initial filtering of dataset
    if args.sex_karyotypes:
        study_sex_karyotypes = set(args.sex_karyotypes.split(','))
        meta = initial_filter(meta,sex_karyotypes=study_sex_karyotypes,cohorts = args.cohorts)
    else:
        meta = initial_filter(meta,cohorts = args.cohorts)

    # Remove samples with irrelevant cancer diagnosis for this study
    sample_size4 = len(meta)

    with open(args.log_file, "a") as f:
        f.write(f"{sample_size4}\t{(sample_size3 - sample_size4)}\t{round(((sample_size3 - sample_size4)/sample_size3),3) * 100}\tRemove samples that are not {args.sex_karyotypes} and not in {args.cohorts} cohort(s).\n")

    
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
        f.write(f"{sample_size5}\t{(sample_size4 - sample_size5)}\t{round(((sample_size4 - sample_size5)/sample_size4),3) * 100}\tExcluded due to relatedness with other individuals.\n")

    # Downsample controls to match on cases ancestry
    min_ancestry_case_control_ratio = math.ceil(get_min_case_control_ratio(meta))
    meta = match_controls_by_ancestry(meta, max_controls_per_case = min_ancestry_case_control_ratio)
    sample_size6 = len(meta)
    with open(args.log_file, "a") as f:
        f.write(f"{sample_size6}\t{(sample_size5 - sample_size6)}\t{round(((sample_size5 - sample_size6)/sample_size5),3) * 100}\tMatched by ancestry allowing 3x controls as cases per ancestry.\n")

    ## Print Summary Statistics
    with open(args.log_file, 'a') as f:
        f.write("===== Summary Report =====\n\n")

        # Total number of individuals
        total_count = len(meta)
        f.write(f"Total individuals: {total_count}\n\n")

        # Count per cohort
        f.write("Counts of cohort per cancer:\n")
        cohort_counts = meta.groupby('cancer')['cohort'].value_counts().unstack(fill_value=0)
        f.write(cohort_counts.to_string())
        f.write("\n\n")

        # Count per cancer
        f.write("Counts by cancer:\n")
        f.write(f"{args.cancer_subtype}\t{(meta['cancer'] != 'control').sum()}\n")
        f.write(f"Controls\t{(meta['cancer'] == 'control').sum()}\n")
        f.write("\n\n")

        # Count per sex_karyotype
        f.write("Counts by sex_karyotype:\n")
        f.write(meta['sex_karyotype'].value_counts().to_string())
        f.write("\n\n")

        # Mean age for case and control
        f.write("Mean age by cancer type:\n")
        mean_case_age = round(meta[meta['cancer'] != "control"]['age'].mean(),2)
        mean_control_age = round(meta[meta['cancer'] == "control"]['age'].mean(),2)
        f.write(f"Mean Case Age\t{mean_case_age}\n")
        f.write(f"Mean Control Age\t{mean_control_age}\n")
        f.write("\n\n")

        # intake_qc_pop count per cancer
        f.write("Counts of intake_qc_pop per cancer:\n")
        case_cohort = meta[meta['cancer'] != "control"]
        control_cohort = meta[meta['cancer'] == "control"]
        case_intake_qc_pop = case_cohort['intake_qc_pop'].value_counts().to_string()
        control_intake_qc_pop = control_cohort['intake_qc_pop'].value_counts().to_string()
        f.write("Cases:\n")
        f.write(case_intake_qc_pop + "\n\n")
        f.write("Controls:\n")
        f.write(control_intake_qc_pop + "\n\n")

        # All cancer types included in this diagnoses are listed below
        f.write(meta['cancer'].value_counts().to_string())
        f.write("\n\n")

        f.write("==========================\n\n")


    meta = meta.merge(pca, left_on='original_id', right_on='#IID', how='left')

    # Write to outfile
    meta.to_csv(args.outfile, sep='\t', index=False, na_rep='NA')

if __name__ == '__main__':
    main()
