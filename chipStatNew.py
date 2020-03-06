#!/usr/bin/env python
import numpy as np
from time import time
from collections import Counter, defaultdict
import sys, os
import resource
import scipy.stats as sc
from hw import HW1

if len(sys.argv) < 3:
    print "Usage: chipStat.py <ped file> <map file>"
    exit()

pedFn = sys.argv[1]
print >>sys.stderr,  pedFn
mapFn = sys.argv[2]
print >>sys.stderr,  mapFn
#sscPosFn = sys.argv[3]
#print >>sys.stderr, sscPosFn
#chunkN = pedFn.split("/")[1].split(".")[0]
#print >>sys.stderr,  chunkN

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

bases=list('ACGT0')

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
famN = 0
CompM = {x:y for x,y in zip(list('ACGT0ID'),list('TGCA0ID'))}
chrM = {'chr'+str(v):i+1 for i,v in enumerate(range(1,23) + ['X'])}
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
print >>sys.stderr, "Read chip.map", time() - st

MAP = [[x[0], x[3]] for x in snpMap]
mapLen = len(MAP)

#GGP = np.array([Counter() for k in range(mapLen)])
#GGC = np.array([Counter() for k in range(mapLen)])

GGP = np.zeros((mapLen,5), dtype=int)
GGC = np.zeros((mapLen,5), dtype=int)

IDX = np.zeros((mapLen), dtype=int)
ALLMap = np.array([defaultdict(int) for k in range(mapLen)])

global start
start = True
positions = defaultdict(list)

class P():
    pass

def processFam(fams):
    for fid, fam in fams.items():
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
            print >>sys.stderr, 'processing family', cur_fam
            if processFam(fams):
                statFam(fams)
            else:
                print >>sys.stderr, "strange family", str(fams)
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

# process the last family
print >>sys.stderr, 'last family', cur_fam
if processFam(fams):
    statFam(fams)
else:
    print sys.stderr, "strange family", str(fams)

statsP = GGP
statsC = GGC

N = parents

for k in range(statsP.shape[0]):
    cnt = sum(statsP[k])
    if  cnt !=  N:
        print >>sys.stderr, k, statsP[k], cnt, "!=", N, parents


#hdb_stat = open('STAT/'+chunkN+'-statNew.npz', 'wb')
hdb_stat = open('chip-statSummary.npz', 'wb')
np.savez(hdb_stat, ALLMap=ALLMap, statsP=statsP, statsC=statsC)
hdb_stat.close()

print >>sys.stderr, "Done"

    
