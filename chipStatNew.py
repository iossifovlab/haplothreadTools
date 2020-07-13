#!/usr/bin/env python
import numpy as np
from time import time
from collections import Counter, defaultdict
import sys, os
import resource
import scipy.stats as sc
from hw import HW1

if len(sys.argv) < 3:
    print("Usage: chipStatNew.py <ped file> <map file>")
    exit()

print(sys.version, 
pedFn = sys.argv[1]
print(pedFn, file=sys.stderr)
mapFn = sys.argv[2]
print(mapFn, file=sys.stderr)

"""
chr  
pos  
refA  
majorA  
minorA  
MM_N  
Mm_N  
mm_N 
00_N  
majorA_N = (2*MM_N + Mm_N) 
minorA_N = (2*mm_N + Mm_N) 
0_N      = (2*OO_N) 
percentCalled = 100*(majorA_N + minorA_N)/(majorA_N + minorA_N + 0_N)  
MAF (minor allele frequency) = minorA_N / (majorA_N + minorA_N) 
HW = hw_test(MM_N,Mm_N,mm_N) 
"""

n=-1
global parents
global children
global GGP
global GGC
parents = 0
children = 0

snpMap = []
persons = {}
nucFams = defaultdict(list)
Gmap={'1':'M', '2':'F', '0':'0'}
S = {'mom':0,'dad':1, 'prb':2,'sib':3}
cur_fam = None
FAM = {}
snpDict = defaultdict()
st = time()
with open(mapFn, 'r') as f:
    n = 0
    for l in f:
        cs = l.strip('\n\r').split("\t")
        ch = 'chr'+cs[0]
        p = cs[3]
        snpMap.append(cs)
        snpDict[(ch,p)] = 1
        n +=1
print("Read chip.map", time() - st, file=sys.stderr)

MAP = [[x[0], x[3]] for x in snpMap]
mapLen = len(MAP)

#
# Since the chip has only two possible bases in each position,
# all possible choices of counts can be represented by 5 element array
# [#(b1 b1), #(b2 b2), #(b1 b2),  #(b2 b1), #(0 0)]
# We do not know in advance which b1 or b2 is the major, 
# so one of #(b1 b2) or #(b2 b1) may be zero.
# In reality only one of (b1 b2) (b2 b1) is present in the ped file,
# so really 4 element array suffies. We still have 5, but the index will 
# start from 1.
#

GGP = np.zeros((mapLen,5), dtype=int)
GGC = np.zeros((mapLen,5), dtype=int)

IDX = np.zeros((mapLen), dtype=int)
ALLMap = np.array([defaultdict(int) for k in range(mapLen)])

global start

class P():
    pass

def processFam(fams):
    for fid, fam in list(fams.items()):
        mom = None
        dade = None
        if len(fam) <3:
            return False
        for v in fam:
            p = P();
            p.mom = v[3]
            if p.mom != '0':
                mom = p.mom
                p.role = 'prb' if v[5] == '2' else 'sib'                
            p.dad = v[2]
            if p.dad != '0':
                dad = p.dad
                p.role = 'prb' if v[5] == '2' else 'sib'
            p.gender = Gmap[v[4]]
            p.personId = v[1]
            p.famId = v[0]
            p.aff = v[5]
            nucFams[fid].append(p)
            persons[p.personId] = p
        if not mom or not dad:
            return False        
        for p in nucFams[fid]:
            if p.personId == mom:
                p.role = 'mom'
            if p.personId == dad:
                p.role = 'dad'
        nucFams[fid] = sorted(nucFams[fid], key=lambda x: S[x.role])
    return True


def CNT(k, A):
    """ For each position k in the chip.map and genotype A=(a,b)
    CNT(k, A) returnes the index in range(5) of the genotype A
    for the position k. The first indext will be 1, for each new genotype
    index is incremented by 1. Indeces are retained in the ALLMap.
    """
    if A in ALLMap[k]:
        return ALLMap[k][A]
    else:
        IDX[k] += 1
        ALLMap[k][A] = IDX[k]
        return ALLMap[k][A]

T={i:[0,0,0,0,0] for i in range(5)}
for i in range(5):
    T[i][i]=1

def statFam(fams):
    global parents
    global children
    global GGP
    global GGC
    members = nucFams[cur_fam]
    for p in members:
        FAM[p.personId]['role']=p.role

    for p in members[:2]:
        parents +=1
        g = FAM[p.personId]['genotype']
        
        GGP += np.array([T[CNT(k, (g[2*k],g[2*k+1]))] for k in range(mapLen)])

    for p in members[2:]:
        children +=1
        g = FAM[p.personId]['genotype']
        GGC += np.array([T[CNT(k, (g[2*k],g[2*k+1]))] for k in range(mapLen)])

n = 0
fams = defaultdict(list)
with open(pedFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        fId, pId, faId, moId, sex, aff = cs[:6]
        #print >>sys.stderr, fId, pId
        if cur_fam and fId != cur_fam:
            #print 'A'
            print('processing family', cur_fam, file=sys.stderr)
            if processFam(fams):
                statFam(fams)
            else:
                print("strange family", str(fams), file=sys.stderr)
            FAM = {}
            fams = defaultdict(list)
        #print 'B'
        cur_fam = fId
        fams[fId].append([fId, pId, faId, moId, sex, aff])
        G = cs[6:]
        FAM[pId] = {'personId':pId,'genotype':np.array(G)}
        n +=1
        #if n == 19:
        #    break

# process last family
print('last family', cur_fam, file=sys.stderr)
if processFam(fams):
    statFam(fams)
else:
    print(sys.stderr, "strange family", str(fams))

statsP = GGP
statsC = GGC

N = parents

for k in range(statsP.shape[0]):
    cnt = sum(statsP[k])
    if  cnt !=  N:
        print(k, statsP[k], cnt, "!=", N, parents, file=sys.stderr)

hdb_stat = open('chip-statSummary.npz', 'wb')
np.savez(hdb_stat, parents=np.array(parents), ALLMap=ALLMap, statsP=statsP, statsC=statsC)
hdb_stat.close()

print("Done", file=sys.stderr)

    
