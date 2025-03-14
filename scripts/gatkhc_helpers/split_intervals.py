#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# The Germline Genomics of Cancer (G2C)
# Copyright (c) 2025-Present, Ryan L. Collins and the Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0

"""
Divide a GATK-style interval_list file into shards, either yielding shards matched
on a fixed target genomic size or balanced for the density of variants in a reference
dataset, like gnomAD
"""


# Import libraries
import argparse
import numpy as np
import pandas as pd
import pybedtools as pbt
from sys import stdin, stdout
from pysam import TabixFile


def split_by_density(int_list, var_beds, vars_per_shard, min_size=5000, verbose=False):
    """
    Splits intervals in int_list into shards where each shard carries no more than
    vars_per_shard total variants present in var_beds
    """

    # Convert intervals to pbt.BedTool and ensure coordinate-sorted
    int_bt = pbt.BedTool('\n'.join(int_list), from_string = True).sort()

    # Open tabix handles to var_beds
    var_tbx = {i : TabixFile(var_beds[i]) for i in range(len(var_beds))}

    # Iterates over each original input interval
    shards = []
    for orig_int in int_bt:
        chrom = orig_int.chrom
        istart = orig_int.start
        iend = orig_int.end

        # Open handles to var beds for positions within this original interval
        ivar_tbxs = {i : var_tbx[i].fetch(chrom, istart, iend) for i, b in var_tbx.items()}

        # Initialize coordinate pointers for this shard
        def _increment_tbx(tbx):
            try:
                line = next(tbx)
                return int(line.split('\t')[1])
                # Once the tabix handle is exhausted, emit an impossibly large coordinate
            except:
                return 10e10
        tbx_pointers = {i : _increment_tbx(t) for i, t in ivar_tbxs.items()}

        # Function to consume current minimum position, increment that pointer, and return updated tbx_pointers
        def _spend_variant(tbx_pointers):
            i = np.argmin(list(tbx_pointers.values()))
            tbx_pointers[i] = _increment_tbx(ivar_tbxs[i])
            return tbx_pointers

        # Count variants until either vars_per_shard have been allocated
        # or until iend is reached, in which case, always emit shard and start new
        s_start = istart
        cur_pos = np.min(list(tbx_pointers.values()) + [iend])
        k = 0
        while cur_pos <= iend:

            # If current position is sufficiently close to end of interval s/t
            # a new shard larger than min_size cannot fit between cur_pos and iend,
            # round cur_pos up to iend and break this loop
            if cur_pos >= iend - min_size:
                cur_pos = iend
                break

            # Find next variant, increment that pointer, and add to tally
            tbx_pointers = _spend_variant(tbx_pointers)
            k += 1
            cur_pos = np.min(list(tbx_pointers.values()) + [iend])

            # Once the number of variants has been tallied, emit current shard
            # and update new shard start/end
            if k >= vars_per_shard:
                shards.append([chrom, s_start, cur_pos])
                if verbose:
                    msg = 'Split interval {}:{:,}-{:,} at {}:{:,}-{:,} ' + \
                          '({:.0f} kb) after {:,} variants\n'
                    stdout.write(msg.format(chrom, istart, iend, chrom, s_start, 
                                            cur_pos, (cur_pos - s_start) / 1000, k))
                cur_pos += 1
                s_start = cur_pos
                cur_pos += 1
                k = 0

        # As the above loop should always break with imperfect sharding, we
        # always need to emit the final shard to avoid right-truncation
        shards.append([chrom, s_start, int(np.max([cur_pos, iend]))])
        if verbose:
            msg = 'Yielding {}:{:,}-{:,} ({:.0f} kb) as the remainder of ' + \
                  'interval {}:{:,}-{:,} after {:,} variants\n'
            stdout.write(msg.format(chrom, s_start, cur_pos, (cur_pos - s_start) / 1000, 
                                    chrom, istart, iend, k))


    return sorted(shards, key = lambda x: int(x[1]))


def split_by_size(int_list, target_size):
    """
    Recursively splits the largest interval in half until mean interval is ≤ --target-size
    """

    int_df = pd.DataFrame(int_list, columns='isize fields'.split())

    avg_size = int_df.isize.mean()
    while avg_size > target_size:

        # Find largest interval
        big_idx = int_df.isize.argmax()
        parent = int_df.fields[big_idx]
        parent[1] = int(parent[1])
        parent[2] = int(parent[2])

        # Divide interval in half
        midpoint = int(np.floor(np.mean([parent[2], parent[1]])))
        children = [(midpoint - parent[1], [*parent[0:2], midpoint, *parent[3:]], ),
                    (parent[2] - midpoint, [parent[0], midpoint+1, *parent[2:]], ),]

        # Drop parent and append children
        int_df = pd.concat([int_df.drop(big_idx),
                            pd.DataFrame(children, columns=int_df.columns)]).\
                    reset_index(drop=True)

        # Recompute mean interval size
        avg_size = int_df.isize.mean()

    # Return list of new intervals sorted by start coordinate
    return sorted(int_df.fields.tolist(), key=lambda x: int(x[1]))


def main():
    """
    Main block
    """
    parser = argparse.ArgumentParser(
             description=__doc__,
             formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-i', '--input-intervals', required=True, 
                        help='Path to input .interval_list', default='stdin')
    parser.add_argument('-t', '--target-size', help='Desired average interval ' +
                        'size. By default, intervals will be divided based on ' +
                        'their genomic size, controlled by this parameter.',
                        type=float, default=10e10)
    parser.add_argument('--n-shards', type=int, help='Total number' +
                        'of desired shards. Only used if --vars-per-shard and ' +
                        '--var-sites are also specified, otherwise will default ' +
                        'to --target-size')
    parser.add_argument('--var-sites', action='append', help='One or more BED ' +
                        'files listing known variant sites in an external ' +
                        'reference dataset, like gnomAD. May be provided ' +
                        'multiple times. If specified with --vars-per-shard, ' +
                        'intervals will be divided into shards balanced by ' +
                        'number of variants per shard.')
    parser.add_argument('--vars-per-shard', type=float, help='Desired number' +
                        'of variants present in --var-sites to allocate per shard. ' +
                        'Required for variant density-based shard balancing ' +
                        'if --var-sites is also optioned.')
    parser.add_argument('--min-interval-size', type=int, default=5000,
                        help='Size of smallest interval shard to be emitted, in ' +
                        'base pairs.')
    parser.add_argument('--gatk-style', action='store_true', default=False,
                        help='Write output as simple GATK-style intervals ' +
                        '[default: write Picard-style intervals with header]')
    parser.add_argument('-o', '--output-intervals', required=True,
                        help='Path to output .interval_list', default='stdin')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print ' +
                        'diagnostic logging to stdout')
    args = parser.parse_args()

    # Determine & report desired splitting mode 
    if args.vars_per_shard is not None \
    and len(args.var_sites) > 0:
        split_mode = 'density'
        msg = 'Splitting intervals by variant density due to the presence of ' + \
              'both --var-sites and --vars-per-shard'
    else:
        split_mode = 'size'
        msg = 'Splitting intervals by size due to the lack of other specific options'
    print(msg + '\n')

    # Open connection to input file
    if args.input_intervals in '- stdin /dev/stdin'.split():
        infile = stdin
    else:
        infile = open(args.input_intervals)

    # Open connection to output file
    if args.output_intervals in '- stdout /dev/stdout'.split():
        outfile = stdout
    else:
        outfile = open(args.output_intervals, 'w')

    # Read input file into memory
    int_list = []
    for line in infile.readlines():

        # Write header from input to output
        if line.startswith('@'):
            if not args.gatk_style:
                outfile.write(line)

        # Store input intervals as tuples of (interval size, list of fields)
        # or as strings for pbt.BedTool, depending on split_mode
        else:
            fields = line.rstrip().split()
            if split_mode == 'density':
                int_list.append('\t'.join(fields[:3]))
            else:
                size = int(fields[2]) - int(fields[1])
                int_list.append((size, fields, ))

    # After reading input intervals, perform splitting based on value of split_mode
    if split_mode == 'density':
        shards = split_by_density(int_list, args.var_sites, args.vars_per_shard, 
                                  args.min_interval_size, args.verbose)
    else:
        shards = split_by_size(int_list, args.target_size)

    # Write sharded intervals to output file
    for vals in shards:
        outline = '{}:{}-{}'.format(*vals[:3])
        if not args.gatk_style:
            outline += '\t' + '\t'.join('+ . intersection ACGTmer'.split())
        outfile.write(outline + '\n')

    # Clear buffer
    outfile.close()


if __name__ == '__main__':
    main()

