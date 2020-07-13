#!/usr/bin/env python
import sys, os
from GenomeAccess import openRef

if len(sys.argv) < 3:
    print("Usage: ref2map.py <map file> <reference file>")
    sys.exit(1)

sep="\t"
mapFn = sys.argv[1]
refFn = sys.argv[2]

GA = openRef(refFn)
chrM = {'chr'+str(v):i+1 for i,v in enumerate(range(1,23))}
chrM['chrX']=23

with open(mapFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split(sep)
        ch = 'chr'+'X' if cs[0] == '23' else cs[0]
        p = cs[3]
        refA = GA.getSequence(ch, int(p),int(p)) if ch in chrM else ''
        print("\t".join(cs + [refA]))

print("Read chip.map", file=sys.stderr)




