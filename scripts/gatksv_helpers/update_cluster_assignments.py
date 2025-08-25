#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Update SV cluster assignments from a second cluster file

Note: this is a helper script intended to be executed within CollapseRedundantSvs.wdl

Functionality for broad, flexible reuse is minimal
"""


import argparse
import networkx as nx
import pybedtools as pbt
from sys import stdout


def add_clusters_to_graph(bedtool, G, record_coords):
    """
    Add nodes and edges from a BED5 cluster assignment file.
    Nodes = record IDs, Edges = cluster membership.
    Disclaimer: this function was written by ChatGPT
    """

    clusters = {}

    for f in bedtool:
        record_coords.setdefault(f.name, (f.chrom, f.start, f.end))
        clusters.setdefault(f[4], []).append(f.name)

    for members in clusters.values():

        # Always add all nodes to the graph
        for rid in members:
            G.add_node(rid)
        
        # Then add edges for all pairs in this cluster
        for i in range(len(members)):
            for j in range(i + 1, len(members)):
                G.add_edge(members[i], members[j])


def main():
    """
    Update cluster assignments by enforcing single-linkage between
    raw and updated BED5 clusterings
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--initial-assignments', required=True, 
                        help='Initial cluster assignments', metavar='BED5')
    parser.add_argument('-u', '--update-assignments', required=True, 
                        help='Cluster assignment updates', metavar='BED5')
    parser.add_argument('-o', '--output-bed', required=True, metavar='BED4',
                        help='Path to output BED4 of collapsed clusters', 
                        default='stdout')
    args = parser.parse_args()

    # Read initial and updated assignments
    raw = pbt.BedTool(args.initial_assignments)
    updates = pbt.BedTool(args.update_assignments)

    # Graph: nodes = record IDs, edges = cluster membership
    G = nx.Graph()

    # Store BED fields for output (record ID → (chrom, start, end))
    record_coords = {}

    # Add clusters from both sources
    add_clusters_to_graph(raw, G, record_coords)
    add_clusters_to_graph(updates, G, record_coords)

    # Assign new cluster IDs based on connected components
    cluster_id = {}
    for cid, component in enumerate(nx.connected_components(G)):
        for rid in component:
            cluster_id[rid] = str(cid)

    # Open connection to output BED4
    if args.output_bed in 'stdout /dev/stdout -'.split():
        fout = stdout
    else:
        fout = open(args.output_bed, 'w')

    # For each cluster, write maximal coordinates and all cluster member IDs
    for cid, component in enumerate(nx.connected_components(G)):
        members = list(component)
        if len(members) < 2:
            continue
        chroms = {record_coords[rid][0] for rid in members}
        if len(chroms) != 1:
            exit(f'Warning: cluster {cid} spans multiple chromosomes: {chroms}\n')
        chrom = list(chroms)[0]
        min_start = min([record_coords[rid][1] for rid in members])
        max_end = max([record_coords[rid][2] for rid in members])
        ids = ','.join(sorted(members))
        fout.write('\t'.join(map(str, [chrom, min_start, max_end, ids])) + '\n')

    # Clear buffer
    if fout is not stdout:
        fout.close()


if __name__ == '__main__':
    main()

