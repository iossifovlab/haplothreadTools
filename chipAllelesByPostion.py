#!/usr/bin/env python
import numpy as np
from time import time
from collections import Counter, defaultdict
import sys, os
import pickle
import resource

if len(sys.argv) < 2:
    print("Usage: chipStat.py <ped file>") 
    exit()


pedFn = sys.argv[1]
print(pedFn, file=sys.stderr)


alleles = None
with open(pedFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        G = cs[6:]
        assert len(G)%2 == 0
        L = len(G)/2
        if alleles:
            assert L == len(alleles)
        else:
            alleles = [defaultdict(int) for p in range(L)]

        for i,g in enumerate(G):
            p = i/2
            alleles[p][g] += 1 

for pd in alleles:
    cs = []
    for a,n in sorted(list(pd.items()),key=lambda x:-x[1]):
        if a != '0':
            cs += [a]
    # if '0' not in pd: cs += ['0','0']
    print("\t".join(cs))
