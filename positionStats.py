#!/usr/bin/env python
import os, sys
from famSummary import *
from collections import defaultdict, Counter
#import pickle
import numpy as np

fams = list(cleanFMS.keys())

if len(sys.argv) < 2:
    print("Usage positionStats.py <chunk #>")
    exit()

chunk = sys.argv[1]
sline = int(chunk)*500
eline = (int(chunk)+1)*500

pedFn =  'SPARK.30K.genotype.release_set.ped'
print(sline, eline, pedFn, "number of parents:", len(fams)*2, file=sys.stderr)

ped = []
n=-1
parentsQ = 0
parentsT = 0
parentsM = 0
GGPQ = defaultdict(list)
GGPT = defaultdict(list)
GGPM = defaultdict(list)
GGCQ = defaultdict(list)
GGCT = defaultdict(list)
GGCM = defaultdict(list)
with open(pedFn, 'r') as f:
    for l in f:
        n +=1
        if n < sline:
            continue
        if n >= eline:
            break
        cs = l.strip('\n\r').split(' ')
        if not cs[0] in fams:
            continue
        pId = cs[1]
        
        member = [[x[0],x[6]] for x in cleanFMS[cs[0]]
                  if x[6] in ["Mother","Father"] and x[0] == pId]
        G = cs[6:]
        if member:
            if cs[0] in quads:
                parentsQ += 1
                {GGPQ[k].append((G[2*k],G[2*k+1])) for k in range(int(len(G)/2))}
            elif cs[0] in trios:
                parentsT += 1
                {GGPT[k].append((G[2*k],G[2*k+1])) for k in range(int(len(G)/2))}
            elif cs[0] in multiplex:
                parentsM += 1
                {GGPM[k].append((G[2*k],G[2*k+1])) for k in range(int(len(G)/2))}                
        else:
            if cs[0] in quads:
                {GGCQ[k].append((G[2*k],G[2*k+1])) for k in range(int(len(G)/2))}
            elif cs[0] in trios:
                {GGCT[k].append((G[2*k],G[2*k+1])) for k in range(int(len(G)/2))}
            elif cs[0] in multiplex:
                {GGCM[k].append((G[2*k],G[2*k+1])) for k in range(int(len(G)/2))} 
            
statsPQ = np.array([Counter(GGPQ[k]) for k in GGPQ])
statsPT = np.array([Counter(GGPT[k]) for k in GGPT])
statsPM = np.array([Counter(GGPM[k]) for k in GGPM])

statsCQ = np.array([Counter(GGCQ[k]) for k in GGCQ])
statsCT = np.array([Counter(GGCT[k]) for k in GGCT])
statsCM = np.array([Counter(GGCM[k]) for k in GGCM])

#fn = open('STATS/posStats-'+chunk+'.pckl', 'wb')
fn = open('STATS/posStats-'+chunk+'.npz', 'wb')
#pickle.dump([statsPQ, statsPT,statsPM,statsCQ,statsCT,statsCM, parentsQ, parentsT, parentsM], fn)
np.savez(fn, SPQ=statsPQ, SPT=statsPT, SPM=statsPM, SCQ=statsCQ, SCT=statsCT, SCM=statsCM, pq=np.array(parentsQ), pt=np.array(parentsT), pm=np.array(parentsM), allow_pickle=True)
fn.close()

print("Done", 'STATS/posStats-'+str(chunk)+'.npz', file=sys.stderr) 
