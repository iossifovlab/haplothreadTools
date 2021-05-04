#!/usr/bin/env python
import os, sys
from collections import defaultdict, Counter
import pickle
import numpy as np

#qFn = 'quadsStruct.txt'
#qFn = 'multiplexStruct.txt'

if len(sys.argv) <3:
    print("Usage: cSummary.py <family structure file> <minFamSize>")
    exit()
    
qFn = sys.argv[1]
minFamSize = int(sys.argv[2])

MZ=[]
with open('MZ.txt', 'r') as f:
    for l in f:
        MZ.append(l.strip('\n\r'))

famStruct = defaultdict(list)
with open(qFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        if minFamSize == 4 and cs[1] in MZ:
            continue
        famStruct[(cs[0],cs[1])] = cs[2]

fams = sorted(set([x[1] for x in famStruct]))

counts = defaultdict()
st = True
posNum = None
for n in range(56):
    with open('STATS/cCnt-'+str(n)+'.npz', 'rb') as f:
        data = np.load(f,encoding='ASCII')
        data.allow_pickle=True
        cnt = data.get('counts')
        for x in cnt:
            if st:
                posNum = float(sum([int(v) for v in list(x[2].values())]))
                st = False
            counts[(x[0],x[1])] = x[2]['0']

fms = [sorted([ [k[1],famStruct[k],k[0]]  for k in famStruct if k[1]==f]) for f in fams]
result = [[ [x[0],x[1],x[2],counts[(x[0],x[2])]] for x in y] for y in fms]
for res in result:
    for r in res:
        r.append(r[-1]/(posNum))

"""
print '\t'.join('famId mom dad aff unaff'.split(' '))
for r in result:
    print '\t'.join(map(str,[r[0][0], r[1][1], r[0][1], r[2][1],r[3][1]))
print '\t'.join('famId role pId #00 fraction'.split(' '))
"""
#"""
fams = defaultdict(list)
for res in result:
    for x in res:
        if x[4] < 0.05:
            fams[x[0]].append((x[1].split('.')[1], x[2]))

clean = [x for x in fams if 'mo' in [p[0] for p in fams[x]] and 'fa' in [p[0] for p in fams[x]]]
for f in sorted(clean):
    if len(fams[f]) >= minFamSize:
        for p in fams[f]:
            print('\t'.join([f,p[1]]))
#"""
