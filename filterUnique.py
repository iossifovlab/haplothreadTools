#!/usr/bin/env python
import os, sys
from collections import defaultdict

#fn = 'statsSummary.txt'
if len(sys.argv) <2:
    print("Usage: filterUnique.py <in file>")
    exit()
    
fn = sys.argv[1]

stats  = defaultdict(int)
with open(fn, 'r') as f:
    f.readline()
    for l in f:
        stats[tuple(l.strip('\n\r').split('\t')[:2])] +=1

with open(fn, 'r') as f:
    cs = f.readline().strip('\n\r').split('\t')
    print('\t'.join(cs))
    for l in f:
        cs = l.strip('\n\r').split('\t')
        if stats[tuple(cs[:2])] == 1:
            print('\t'.join(cs))


