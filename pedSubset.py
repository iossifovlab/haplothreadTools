#!/usr/bin/env python
import sys,os
import numpy as np

if len(sys.argv) < 4:
    print("Usage: pedSubset.py <pos indices> <input ped> <out ped>")
    exit()


posFn = sys.argv[1]
inPed = sys.argv[2]
outPed = sys.argv[3]

pos = []
with open(posFn, 'r') as f:
    for l in f:
        pos.append(l.strip('\n\r'))

pos = list(map(int,pos))
P = np.array(sorted(list(range(6)) + [6+2*n for n in pos] + [6+2*n+1 for n in pos]))

print(P.shape, file=sys.stderr)


with open(outPed, 'w') as out:
    with open(inPed, 'r') as f:
        for l in f:
            cs = np.array(l.strip('\n\r').split(' '))
            cs = cs[P]
            out.write(' '.join(cs)+'\n')


    
    
