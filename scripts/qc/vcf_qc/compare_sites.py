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


# Declare common constants
# Note that short variant subclasses are lowercase and structural variant
# subclasses are uppercase to distinguish them by size
var_classes = 'snv indel sv'.split()
var_subclasses = {'snv' : 'ti tv'.split(),
                  'indel' : 'ins del'.split(),
                  'sv' : 'DEL DUP CNV INS INV CPX CTX'.split()}
vsc_matches = {'del' : 'del DEL CNV'.split(),
               'ins' : 'ins INS DUP CNV'.split(),
               'DEL' : 'del DEL CNV'.split(),
               'DUP' : 'ins INS DUP CNV'.split(),
               'CNV' : 'ins del DEL DUP CNV'.split(),
               'INS' : 'ins INS DUP CNV'.split(),
               'INV' : 'INV CPX'.split(),
               'CPX' : 'INV CPX'.split(),
               'CTX' : 'CTX'}


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


def find_overlap_hits(hits_g, a_bt, b_bt, genome, ro=0.1, min_size=10, 
                      max_size_diff=3.0, max_bkpt_rel=0.5, pad=1):
    """
    Search for matching variants between two pbt.BedTool objects based on overlap
    """

    # Filter each pbt.BedTool to variants larger than min_size
    # Note that here we want to use the variant size encoded in the site metrics
    # file, not the pbt.Interval.size (to handle cases of 1bp intervals for INS/CTX)
    big_nids = [nid for nid in hits_g.nodes() 
                if hits_g.nodes[nid].get('size', 0) >= min_size]
    a_bt = a_bt.filter(lambda x: 'a_' + x[3] in big_nids).\
                cut(range(6)).slop(b=pad, g=genome)
    b_bt = b_bt.filter(lambda x: 'b_' + x[3] in big_nids).\
                cut(range(6)).slop(b=pad, g=genome)

    # Find overlaps, filter, and add qualifying edges to graph
    for ovr in a_bt.intersect(b_bt, wao=True, f=ro, r=True):

        # Check for any hit
        if ovr[6] == '.':
            continue

        # Get left variant info
        nid_a = 'a_' + ovr[3]
        start_a = hits_g.nodes[nid_a].get('pos')
        end_a = hits_g.nodes[nid_a].get('end')
        size_a = hits_g.nodes[nid_a].get('size', 0)
        vsc_a = hits_g.nodes[nid_a].get('vsc')
        af_a = hits_g.nodes[nid_a].get('af')

        # Get right variant info
        nid_b = 'b_' + ovr[9]
        start_b = int(ovr[7])
        end_b = int(ovr[8])
        size_b = hits_g.nodes[nid_b].get('size', 0)
        vsc_b = hits_g.nodes[nid_b].get('vsc')
        af_b = hits_g.nodes[nid_b].get('af')

        # Check for matching variant types
        if vsc_b not in vsc_matches[vsc_a]:
            continue

        # Check for sufficiently similar sizes
        smaller, larger = sorted([size_a, size_b])
        if larger / smaller > max_size_diff:
            continue

        # Check for reasonably close left OR right breakpoints
        left_dist = np.abs(start_a - start_b)
        right_dist = np.abs(end_a - end_b)
        if min([left_dist, right_dist]) / larger > max_bkpt_rel:
            continue

        # Compute Euclidean distance between variants and update hits graph
        e_d = node_distance([start_a, start_b], [end_a, end_b], [af_a, af_b])
        hits_g.add_edge(nid_a, nid_b, dist=e_d)

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

    lbp_d = np.abs(starts[0] - starts[1]) / bp_span

    rbp_d = np.abs(ends[0] - ends[1]) / bp_span

    af_d = np.abs(afs[0] - afs[1])

    return np.sqrt(((1 - ovr_jac) ** 2) + (lbp_d ** 2) + (rbp_d ** 2) + (af_d ** 2))


def get_distances(hits_g, node_id):
    """
    Return a dict of {match: distance}, ordered from closest to farthest
    """

    dists = {}
    for edge in hits_g.edges(node_id):
        other_id = list(set(edge).difference(set([node_id])))[0]
        dists[other_id] = hits_g.edges[edge].get('dist')
    
    return dict(sorted(dists.items(), key=lambda x: x[1]))


def prune_hits(hits_g):
    """
    Prune edges in a graph of overlapping variants until no node has < 2 edges
    """

    node_ids = hits_g.nodes()

    n_edges = {nid : len(hits_g.edges(nid)) for nid in node_ids}

    # Loop and prune edges until all nodes are first-degree
    while max(n_edges.values()) > 1:

        # Select a node of degree > 1 to prune, starting with nodes with the fewest edges
        min_multi_e = np.min([n for n in n_edges.values() if n > 1])
        nid = [nid for nid, n in n_edges.items() if n == min_multi_e][0]

        # Get distance to all candidate matches
        cand_d = get_distances(hits_g, nid)
        cand_ids = list(cand_d.keys())

        # Removing candidates if they are a mutually optimal match to a different variant
        # This ensures variants are always paired with their joint optimal match
        for cid in cand_ids:
            cid_best = list(get_distances(hits_g, cid).keys())[0]
            cid_best_best = list(get_distances(hits_g, cid_best).keys())[0]
            if cid == cid_best_best \
            and cid_best != nid:
                cand_d.pop(cid)

        # Remove all edges other than final candidate, if any remain
        if len(cand_d) == 0:
            final_best_cid = ''
        else:
            final_best_cid = list(cand_d.keys())[0]
        for cid in cand_ids:
            if cid == final_best_cid:
                continue
            hits_g.remove_edge(nid, cid)

        # Update n_edges
        n_edges[nid] = len(hits_g.edges(nid))

    return hits_g


def format_output_bed(hits_g, target_prefix, ref_prefix, genome):
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
    return pbt.BedTool('\n'.join(bt_strs), from_string=True).sort(g=genome)


def compress_overlap_distribs(bt, mode='size', max_size=1000000, 
                              af_d_breaks=[0.01, 0.1, 0.5, 1]):
    """
    Compress variant overlap information contained in a pbt.BedTool 
    by variant class, subclass, and either size or AF bins
    
    This will produce compressed distributions formatted nearly identically to 
    those produced by clean_site_metrics.py for downstream compatability
    """

    # Convert bedtool to dataframe and annotate difference in AF
    df = bt.to_dataframe(header=0, disable_auto_names=True)
    df['af_d'] = np.abs(df.af - df.match_af)

    # Bin absolute difference in allele frequencies
    # The first bin will always correspond to variants with no match
    af_d_le = np.append(-1, np.array(af_d_breaks))
    df['af_d_bin'] = df.af_d.apply(lambda x: np.argmax(x < af_d_le))
    af_d_labels = ['no_match'] + ['af_d_le{:.2e}'.format(float(x)) for x in af_d_breaks]

    # Build pd.DataFrame to collect results
    res_cols = ['#class', 'subclass', mode] + af_d_labels
    res = pd.DataFrame(columns=res_cols)

    # Remaining behavior depends on value of `mode`
    if mode == 'size':
        
        # Assign each variant to a size bin
        size_ge = np.array([0] + [10 ** k for k in 
                                  range(int(np.floor(np.log10(max_size)) + 1))])
        def _assign_size_bin(x):
            return np.argmin(np.append(x >= size_ge, False)) - 1
        df['size_bin'] = df['size'].apply(_assign_size_bin)

        # Iterate over variant classes, subclasses, and size bins
        for vc in var_classes:
            for vsc in var_subclasses[vc]:
                for size_bin in range(len(size_ge)):

                    # Collect variants corresponding to this tranche
                    sub_df = df[(df.vc == vc) & \
                                (df.vsc == vsc) & \
                                (df.size_bin == size_bin)]
                    if len(sub_df) == 0:
                        continue

                    # Count variants in each AF concordance bin
                    out_vals = [vc, vsc, 'ge{}bp'.format(size_ge[size_bin])]
                    for d_i in range(len(af_d_le)):
                        out_vals.append((sub_df.af_d_bin == d_i).sum())
                    res.loc[len(res)] = out_vals

    elif mode == 'af':

        # Assign each variant to an AF bin
        min_af = np.nanmin(df.af)
        af_le = np.array([10 ** k for k in range(int(np.ceil(np.log10(min_af))), 1, 1)])
        df['af_bin'] = df.af.apply(lambda x: np.argmax(x <= af_le))

        # Iterate over variant classes, subclasses, and AF bins
        for vc in var_classes:
            for vsc in var_subclasses[vc]:
                for af_bin in range(len(af_le)):

                    # Collect variants corresponding to this tranche
                    sub_df = df[(df.vc == vc) & \
                                (df.vsc == vsc) & \
                                (df.af_bin == af_bin)]
                    if len(sub_df) == 0:
                        continue

                    # Count variants in each AF concordance bin
                    out_vals = [vc, vsc, 'le{:.2e}'.format(af_le[af_bin])]
                    for d_i in range(len(af_d_le)):
                        out_vals.append((sub_df.af_d_bin == d_i).sum())
                    res.loc[len(res)] = out_vals

    return res


def write_outputs(hits_g, out_prefix, query_prefix, ref_prefix, genome,
                  common_af=None, gzip=False):
    """
    Format and write output files
    """

    # Prep header for sites BED output
    bed_header_fields = '#chrom pos end vid vc vsc size af match_vid match_af dist'
    bed_header = '\t'.join(bed_header_fields.split())

    # BED of outer join
    bt = format_output_bed(hits_g, query_prefix, ref_prefix, genome)
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
    size_d_out = out_prefix + '.concordance_by_size.tsv'
    if gzip:
        size_d_out += '.gz'
    size_d.to_csv(size_d_out, sep='\t', index=False)

    # Compress outer join by variant class, subclass, and AF
    af_d = compress_overlap_distribs(bt, mode='af')
    af_d_out = out_prefix + '.concordance_by_af.tsv'
    if gzip:
        af_d_out += '.gz'
    af_d.to_csv(af_d_out, sep='\t', index=False)


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
    parser.add_argument('-g', '--genome', metavar='file', required=True,
                        help='BEDTools-style .genome file. Required.')
    parser.add_argument('-o', '--output-prefix', help='Prefix for all output files',
                        metavar='path', default='./compare_sites')
    parser.add_argument('-m', '--mode', choices='exact overlap both'.split(),
                        default='exact', help='Criteria to impose when searching ' +
                        'for overlapping variants. [default: exact matches]')
    parser.add_argument('-r', '--min-reciprocal-overlap', type=float, default=0.1,
                        metavar='float', help='Minimum reciprocal overlap to ' +
                        'permit for --mode "overlap" or "both" [default: 0.1]')
    parser.add_argument('-s', '--min-overlap-var-size', type=int, default=10,
                        metavar='int', help='Minimum variant size to ' +
                        'consider for --mode "overlap" or "both" [default: 10]')
    parser.add_argument('--max-overlap-size-diff', type=float, default=3.0,
                        metavar='float', help='Maximum difference in variant sizes ' +
                        'to tolerate for --mode "overlap" or "both" [default: 3.0]')
    parser.add_argument('--max-overlap-bkpt-rel-dist', type=float, default=0.5,
                        metavar='float', help='Maximum distance between either ' +
                        'left or right breakpoints to tolerate for --mode ' +
                        '"overlap" or "both". Specified as a fraction of larger ' +
                        'variant total size. [default: 0.5]')
    parser.add_argument('--overlap-pad', type=int, default=1, metavar='int', 
                        help='Distance to pad each variant for --mode "overlap" ' +
                        'or "both". Does not alter variant coordinates or size ' +
                        'estimates; helps for 0bp/1bp SV intervals like insertions ' +
                        '[default: 1]')
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
    if args.mode in 'exact both'.split():
        hits = find_exact_hits(hits)

    # Look for matches based on overlap, if optioned
    if args.mode in 'overlap both'.split():
        hits = find_overlap_hits(hits, a_bt, b_bt, args.genome,
                                 args.min_reciprocal_overlap, 
                                 args.min_overlap_var_size, 
                                 args.max_overlap_size_diff,
                                 args.max_overlap_bkpt_rel_dist, 
                                 args.overlap_pad)

    # Process candidate matches to generate 1:1 mapping of hits
    hits = prune_hits(hits)

    # Write output files for left outer join  (LoJ; all A sites, with matching B info)
    write_outputs(hits, args.output_prefix + '.loj', 'a_', 'b_', args.genome,
                  args.common_af, args.gzip)

    # Write output files for right outer join  (RoJ; all B sites, with matching A info)
    if not args.no_reverse:
        write_outputs(hits, args.output_prefix + '.roj', 'b_', 'a_', 
                      args.genome, args.common_af, args.gzip)


if __name__ == '__main__':
    main()

