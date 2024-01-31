import math, sys, os
from collections import deque
import numpy as np

if len(sys.argv < 4):
  print("Script requires 3 arguments: input vcf, buffer length, and maximum buffered span")

in_vcf = sys.argv[1]
buffer = int(sys.argv[2])
max_span = int(sys.argv[3])

xrms = np.loadtxt(in_vcf, usecols=(0), dtype='U')
pos = np.loadtxt(in_vcf, usecols=(1), dtype=np.int32)

if len(pos) == 0:
  print("Found empty vcf.")
  os.rename(in_vcf, f'{in_vcf}.subshard_0.vcf')
  sys.exit(0)
    

# unfortunately numpy sorts before returning the unique elements,
#  so we have to unsort (i.e. re-sort the indices into order)
Xorder = list(xrms[sorted(np.unique(xrms, return_index=True)[1])])

L, U = np.clip(pos - buffer, 1, None), pos + buffer             # Lower, Upper buffers around each variant
I = np.nonzero((L[1:] > U[:-1]) | (xrms[1:] != xrms[:-1]))[0]   # where buffer boundaries _dont_ overlap
iA, iB = np.concatenate([[0], I+1]), np.concatenate([I, [-1]])  # lower, upper boundary buffer indices
X, A, B = xrms[iA], L[iA], U[iB]                                # bounding buffers
S = (B-A).sum()                                                 # total buffer span
T = math.ceil(S/math.ceil(S/max_span))                          # target fragment length

print(S, 'sized buffer to', T, 'length fragments')

# Yay borel sets
# deques pop from the end, so we have to load in backwards
borel = deque(list(zip(X, A, B))[::-1])
stop = deque()
s = 0

while len(borel) > 0:
    x, a, b = borel.pop() # deque supports O(1) resizing
    if b-a > (r := T-S): # split a buffer if it overflows a fragment
        m = a + r
        borel.append((x,m,b)) # another reason for the deque
        b = m
    s += b-a
    if s >= T:
        stop.append((x, b-1)) # accumulating fragment stopping points
        s = 0
stop.append((x, b))

print(f'Writing {len(stop)} fragments')
# And now the file work, which is always chunky
with open(in_vcf) as inp:
    header = ''
    for line in inp: # this is gonna bleed onto the first variant line
        if not line.startswith('#'): break
        header += line
    
    i = 0
    out = open(f'{in_vcf}.subshard_{i}.vcf', 'w')
    out.write(header)
    xf, pf = stop[i]
    xfi = Xorder.index(xf)
    out.write(line) # first variant is assumably in first fragment
    
    xrml = '-1'
    xrmi = '-1'
    
    for line in inp:
        xrm, pos, *_ = line.split('\t')
        
        if xrmi != xrm:
            xrml = xrm
            xrmi = Xorder.index(xrm)
        
        if (xrmi > xfi) or (xrmi == xfi and int(pos) > pf):
            out.close()
            
            i += 1
            out = open(f'{in_vcf}.subshard_{i}.vcf', 'w')
            out.write(header)
            xf, pf = stop[i]
            xfi = Xorder.index(xf)
        out.write(line)
        
    out.close()
print('Done!')
