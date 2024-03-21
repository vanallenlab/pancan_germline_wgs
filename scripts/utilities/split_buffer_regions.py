#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#
# Copyright (c) 2024-Present, Kyler K. A. Anderson and the Dana-Farber Cancer Institute
# Contact: Kyler Anderson <Kyler_Anderson@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0
#
# This script takes a vcf, a region span, and a variant buffer as input.
# It measures the (discontinuous) buffer region around the variant positions,
#   then divides it evenly into pieces less than (or rarely equal to) the provided span.
# Specifically, it outputs regions that will divide the input vcf into
#   these evenly sized pieces when used with bcftools +scatter -S;
#   NOT necessarily (and frankly not probably) the pieces themselves. More at the bottom.

import math, sys, os
from collections import deque
import numpy as np

if len(sys.argv) < 4:
  print("Script requires 3 arguments: input vcf, maximum region span, and buffer size.")
  sys.exit(1)

in_vcf = sys.argv[1]
if not in_vcf.endswith('.vcf'):
  print("Input vcf must be in vcf format! Attempting with", os.path.basename(in_vcf))
out_prefix = os.path.splitext(in_vcf)[0]

max_span = int(sys.argv[2])
buffer = int(sys.argv[3])

xrms = np.loadtxt(in_vcf, usecols=(0), dtype='U')
pos = np.loadtxt(in_vcf, usecols=(1), dtype=np.int32)

# unfortunately numpy sorts before returning the unique elements,
#  so we have to unsort (i.e. re-sort the indices into order)
# Xorder = list(xrms[sorted(np.unique(xrms, return_index=True)[1])])

L, U = np.clip(pos - buffer, 1, None), pos + buffer             # Lower, Upper buffers around each variant
I = np.nonzero((L[1:] > U[:-1]) | (xrms[1:] != xrms[:-1]))[0]   # where buffer boundaries _dont_ overlap
iA, iB = np.concatenate([[0], I+1]), np.concatenate([I, [-1]])  # lower, upper merged buffer indices
X, A, B = xrms[iA], L[iA], U[iB]                                # merged buffers (i.e. subregions)
S = (B-A).sum()                                                 # total subregion span
N = math.ceil(S/max_span)                                       # number of regions
T = math.ceil(S/N)                                              # target region length (for equal distribution)

# N*T >= S, so measuring out T will naturally produce N shards.
#  Therefore we advance processively on T, rather than prescriptively on N.
#  Believe it or not, the logic ends up simpler (specifically because of subregions bigger than T).
print(S, 'sized buffer to', T, 'length fragments')

queue = deque(zip(X, A, B))
queue.reverse()             # deques pop from the end, and this is free I think
s, i, N = 0, 0, True        # honest coincidence
with open(f"{out_prefix}.scatter_regions.txt", 'w') as out:
  while len(queue) > 0:
    x, a, b = queue.pop()   # deque supports O(1) resizing
    if N: xl, al = x, a     # if the region is new we need to start a new record
    elif x != xl:           # if the region didn't fill and the chromosome changes...  
      out.write(f'{xl}:{al}-{bl}\t{i}\n') # we need to record what we had...
      xl, al = x, a                       # and start a new record
    if b-a > (r := T-s):    # split a subregion if it overflows a region, r is what will fit
      m = a + r             # mark the new end point to fill the region; half open intervals so m is end and next start
      queue.append((x,m,b)) # the chromosome won't change next pass, but the region will. ...
      b = m                 #   while it's unlikely for this to not happen, it's not impossible
    s += b-a                # contribute subregion measure, half open intervals -> no +1
    if (N := s >= T):       # N signals a new region/this region filled
      out.write(f'{x}:{al}-{b}\t{i}\n')   # x == xl here, so al and b are compatible
      i += 1                # increment region count
      s = 0                 # reset region measure
    bl = b                  # remember where we left off in case the chromosome changes
  if not N: out.write(f'{x}:{al}-{b}\t{i}\n') # If the region didnt fill and there's no next subregion, we need to record what we have
print('Done!')

# Since the output region file is intended only to slice the input vcf, we're not concerned about there being anything between the subregions.
# So rather than record every subregion individually, we "fill in" the gaps between the subregions (merged buffers).
#   xl and al keep track of the beginning of the filled in regions. xl and bl are in case the chromosome changes, which al naturally tracks with.
# The resulting region file will therefore have individual entries that are larger than region_span;
#   when used to slice the input vcf then take +-variant_buffer around each position, the resulting merges will each be the desired size.
