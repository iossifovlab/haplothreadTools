#!/usr/bin/env python
import sys,os

if len(sys.argv) < 4:
    print("Usage: posIndices.py <filteredPos file> <SPARK map in> <SPARK map out> <posIds file>")
    exit()

posFn = sys.argv[1]
mapFn = sys.argv[2]
outMapFn = sys.argv[3]
outIdsFn = sys.argv[4]

pos = {}
with open(posFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        pos[(cs[0],cs[1])] = 1

MAP = []
idx = []
n = -1
with open(mapFn, 'r') as f:
    for l in f:
        n += 1
        cs = l.strip('\n\r').split('\t')
        if (cs[0],cs[3]) in pos:
            MAP.append(l)
            idx.append(n)

with open(outMapFn,'w') as out:
    for l in MAP:
        out.write(l)

with open(outIdsFn, 'w') as f: 
    for id in idx:
        f.write(str(id) + '\n')



