#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Compare two variant site summary .BED files
"""


import argparse
import gzip
import networkx as nx
import numpy as np
import pandas as pd
import pybedtools as pbt
from Bio import bgzf
from collections import Counter
from os import remove
from re import sub


def populate_nodes(hits_g, bt, prefix=''):
    """
    Add all variants present in a pbt.Bedtool (bt) as nodes in a nx.Graph (hits_g)
    """

    for vt in bt:
        vid = prefix + vt.name
        hits_g.add_node(vid, chrom=vt.chrom, pos=int(vt.start), end=int(vt.end), 
                        vc=vt.fields[4], vsc=vt.fields[5], size=int(vt.fields[6]),
                        af=float(vt.fields[8]))

    return hits_g


def find_exact_hits(hits_g):
    """
    Search for exact variant matches of variant ID between two pbt.BedTool 
    objects and updates a graph of hits

    This expects identical variants to have identical IDs, which will be the 
    case if the input BED files were generated by clean_site_metrics.py

    Exact matches are always treated as having distance = 0 irrespective of their AFs
    """

    n_k = Counter([sub('^a_|^b_', '', nid) for nid in hits_g.nodes()])

    match_ids = [nid for nid, k in n_k.items() if k > 1]

    for base_id in match_ids:
        v1 = 'a_' + base_id
        v2 = 'b_' + base_id
        hits_g.add_edge(v1, v2, dist=0.0)

    return hits_g


def node_distance(starts, ends, afs):
    """
    Compute the distance between two variants

    Distance is defined as Euclidean/linear according to:
    A. 1 - footprint jaccard index (bp overlapping / total bp spanned by both variants)
    B. normalized left breakpoint distance (distance between starts / total bp spanned by both variants)
    C. normalized right breakpoint distance (distance between ends / total bp spanned by both variants)
    D. allele frequency difference
    """

    bp_span = np.abs(np.max(ends) - np.min(starts))
    bp_ovr = np.abs(np.min(ends) - np.max(starts))
    ovr_jac = bp_ovr / bp_span

    lbp_d = (np.max(starts) - np.min(starts)) / bp_span

    rbp_d = (np.max(ends) - np.min(ends)) / bp_span

    af_d = np.max(afs) - np.min(afs)

    return np.sqrt((ovr_jac ** 2) + (lbp_d ** 2) + (rbp_d ** 2) + (af_d ** 2))


def prune_hits(hits_g):
    """
    Prune edges in a graph of overlapping variants until no node has < 2 edges
    """

    node_ids = hits_g.nodes()

    n_edges = {nid : len(hits_g.edges(nid)) for nid in node_ids}

    while max(n_edges.values()) > 1:
        # TODO: implement this
        import pdb; pdb.set_trace()

    return hits_g


def format_output_bed(hits_g, target_prefix, ref_prefix):
    """
    Format a pbt.BedTool summarizing overlap results
    """

    node_ids = list(hits_g.nodes())

    # Simplified output line:
    # chrom, start, end, id, vc, vsc, size, af, match_id, match_af, match_dist
    bed_line_fmt = '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2e}\t{}\t{:.2e}\t{:.2}'

    # Gather BED 
    bt_strs = []
    for nid in node_ids:
        # Skip nodes in B file
        if nid.startswith(ref_prefix):
            continue

        # Revert variant ID back to original format
        vid = sub('^' + target_prefix, '', nid)
        ninfo = hits_g.nodes[nid]

        # Gather matching variant information
        edges = list(hits_g.edges(nid))
        if len(edges) == 0:
            match_id = '.'
            match_af = np.nan
            match_dist = np.nan
        elif len(edges) == 1:
            match_nid = list(set(edges[0]).difference(set([nid])))[0]
            match_id = sub('^' + ref_prefix, '', match_nid)
            match_af = hits_g.nodes[match_nid].get('af')
            match_dist = hits_g[nid][match_nid]['dist']
        else:
            msg = 'Node {} from has more than one edge after pruning. ' + \
                  'This indicates a bug that needs to be fixed. Exiting.'
            exit(msg.format(nid))

        # Format output line and add to bt collector
        bed_line_vals = [ninfo[k] for k in 'chrom pos end vc vsc size af'.split()]
        bed_line_vals.insert(3, vid)
        bed_line_vals += [match_id, match_af, match_dist]
        bt_strs.append(bed_line_fmt.format(*bed_line_vals))

    # Return sorted bedtool
    return pbt.BedTool('\n'.join(bt_strs), from_string=True).sort()


def compress_overlap_distribs(bt, mode='size', max_size=1000000):
    """
    Compress variant overlap information contained in a pbt.BedTool 
    by variant class, subclass, and either size or AF bins
    
    This will produce compressed distributions formatted nearly identically to 
    those produced by clean_site_metrics.py for downstream compatability
    """

    # Convert bedtool to dataframe and annotate difference in AF
    df = bt.to_dataframe(header=0, disable_auto_names=True)
    df['af_d'] = np.abs(df.af - df.match_af)

    # 

    # Behavior depends on value of `mode`
    if mode == 'size':
        size_ge = [0] + [10 ** int(k) for k in range(floor(np.log10(max_size)))]

    # Get AF parameters
    min_af = np.nanmin(df.af)
    af_breaks = range(np.ceil(np.log10(min_af)))

    import pdb; pdb.set_trace()


def write_outputs(hits_g, out_prefix, query_prefix, ref_prefix, common_af=None, gzip=False):
    """
    Format and write output files
    """

    # Prep header for sites BED output
    bed_header_fields = '#chrom pos end vid vc vsc size af match_vid match_af dist'
    bed_header = '\t'.join(bed_header_fields.split())

    # BED of outer join
    bt = format_output_bed(hits_g, query_prefix, ref_prefix)
    bed_out = out_prefix + '.sites.bed'
    bt = bt.saveas(bed_out, trackline=bed_header)
    if gzip:
        bt = bt.tabix(force=True)
        remove(bed_out)

    # Subset outer join to common sites and write common subset to file, if optioned
    if common_af is not None:
        common_bed_out = out_prefix + '.sites.common.bed'
        common_bt = bt.filter(lambda f: float(f[7]) >= common_af).saveas()
        common_bt = common_bt.saveas(common_bed_out, trackline=bed_header)
        if gzip:
            common_bt = common_bt.tabix(force=True)
            remove(common_bed_out)

    # Compress outer join by variant class, subclass, and size
    size_d = compress_overlap_distribs(bt, mode='size')
    # TODO: write this result to file

    # # Compress outer join by variant class, subclass, and AF
    # af_d = compress_overlap_distribs(bt, mode='af')
    # # TODO: write this result to file


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-a', '--sites-a', metavar='file', required=True,
                        help='First sites.bed file. Must be sorted. Required.')
    parser.add_argument('-b', '--sites-b', metavar='file', required=True,
                        help='Second sites.bed file. Must be sorted. Required.')
    parser.add_argument('-o', '--output-prefix', help='Prefix for all output files',
                        metavar='path', default='./compare_sites')
    parser.add_argument('-m', '--comparison-mode', choices='exact overlap both'.split(),
                        default='exact', help='Criteria to impose when searching ' +
                        'for overlapping variants. [default: exact matches]')
    parser.add_argument('-r', '--min-reciprocal-overlap', type=float, default=0.1,
                        metavar='float', help='Minimum reciprocal overlap to ' +
                        'permit for --comparison-mode "overlap" or "both" ' +
                        '[default: 0.1]')
    parser.add_argument('-s', '--min-overlap-var-size', type=int, default=10,
                        metavar='int', help='Minimum variant size to ' +
                        'consider for --comparison-mode "overlap" or "both" ' +
                        '[default: 10]')
    parser.add_argument('--common-af', type=float, help='AF cutoff for common ' +
                        'variants. If provided, will generate separate sites.bed ' +
                        'output files restricted to common variants. [default: ' +
                        'do not generate common variant-only .bed outputs]',
                        metavar='[float]')
    parser.add_argument('--no-reverse', action='store_true', help='Do not perform ' +
                        'the complementary analysis when using --sites-b as the ' +
                        'query dataset and --sites-a as the reference [default: ' +
                        'perform analyses in both directions]')
    parser.add_argument('-z', '--gzip', action='store_true', help='Compress ' +
                        'output files with gzip/bgzip [default: write ' + 
                        'uncompressed .tsv/.bed]')
    args = parser.parse_args()

    # Open connections to BED files as pbt.BedTool
    a_bt = pbt.BedTool(args.sites_a)
    b_bt = pbt.BedTool(args.sites_b)

    # Record candidate overlaps as an nx.Graph
    hits = nx.Graph()
    hits = populate_nodes(hits, a_bt, prefix='a_')
    hits = populate_nodes(hits, b_bt, prefix='b_')

    # Look for exact matches, if optioned
    if args.comparison_mode in 'exact both'.split():
        hits = find_exact_hits(hits)

    # Look for matches based on overlap, if optioned
    # TODO: implement this

    # Process candidate matches to generate 1:1 mapping of hits
    hits = prune_hits(hits)

    # Write output files for left outer join  (LoJ; all A sites, with matching B info)
    write_outputs(hits, args.output_prefix + '.loj', 'a_', 'b_', args.common_af, args.gzip)

    # Write output files for right outer join  (RoJ; all B sites, with matching A info)
    if not args.no_reverse:
        write_outputs(hits, args.output_prefix + '.roj', 'b_', 'a_', 
                      args.common_af, args.gzip)


if __name__ == '__main__':
    main()

