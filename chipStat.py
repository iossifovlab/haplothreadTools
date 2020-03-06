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
chunkN = pedFn.split("/")[1].split(".")[0]
print >>sys.stderr,  chunkN

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

"""
sscPos = {}
with open(sscPosFn, 'r') as f:
    for l in f:
        if l.startswith('C'):
            continue
        cs = l.strip('\n\r').split('\t')
        if not (cs[0],cs[1]) in snpDict:
            continue
        cnt = sorted([[bases[i-3], int(cs[i])] for i in range(3,7) if int(cs[i]) > 0], key=lambda x: x[1], reverse=True )
        sscPos[(cs[0],cs[1])] = ''.join([x[0] for x in cnt])

print >>sys.stderr, "Read ssc positions", time() - st
"""

MAP = [[x[0], x[3]] for x in snpMap]
mapLen = len(MAP)
GGP = np.array([Counter() for k in range(mapLen)])
GGC = np.array([Counter() for k in range(mapLen)])

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
        GGP += np.array([Counter([(FAM[p.personId]['genotype'][2*k],
                                   FAM[p.personId]['genotype'][2*k+1])])
                 for k in range(mapLen)])

    for p in members[2:]:
        children +=1
        GGC += np.array([Counter([(FAM[p.personId]['genotype'][2*k],
                                   FAM[p.personId]['genotype'][2*k+1])])
                 for k in range(mapLen)])

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
            if n % 100 == 0:
                print >>sys.stderr, "memory", resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, "time", time() -st , "sec"
                print >>sys.stderr, "time", time() -st , "sec"
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

def verify(sP, k):
    basesP = [item for sublist in sP for item in sublist]
    c = {x:sum([sP[y]*2 if (x,x) == y else sP[y]
                for y in sP if x in y]) for x in basesP}
    return sum(c.values())

def alleles(sP, sC, k):
    basesP = [item for sublist in sP for item in sublist]
    c = {x:sum([sP[y]*2 if (x,x) == y else sP[y]
                for y in sP if x in y]) for x in basesP}
    d = sorted(c.items(), key=lambda x: x[1], reverse=True)
    d = [x for x in d if x[0] != '0']
    if len(d) ==0:
        return ['0', '0', N*2, 0, 0, N*2, 0, 0, N*4, 0,0,0]

    majorA = d[0][0]
    
    if len(d)>1:
        minorA = d[1][0]
    else:
        minorA = majorA
        # 
        # look if  children have another alleles
        basesC = [item for sublist in sC for item in sublist]        
        c = {x:sum([sC[y]*2 if (x,x) == y else sC[y] for y in sC if x in y])
             for x in basesC}
        d = sorted(c.items(), key=lambda x: x[1], reverse=True)
        #print >>sys.stderr, '\t'.join(map(str,[MAP[k][0],MAP[k][1],
        #                                 'parents', sP, 'children', sC, d]))
        d = [x for x in d if x[0] != '0']
        
        if len(d) > 1:
            minorA = d[1][0]
        if len(d) > 2:
            print >>sys.stderr, "Strange", d
    
    keys = [[(majorA, majorA)],
            [(minorA, minorA)],
            [(majorA, minorA), (minorA, majorA)],
            [('0','0')]]

    MM_N, mm_N, Mm_N, OO_N = [sum([sP[i] for i in k]) for k in keys]
    #print >>sys.stderr, majorA, minorA, MM_N, mm_N, Mm_N, OO_N
    if minorA == majorA:
        mm_N = 0
        Mm_N = 0
    majorA_N = 2*MM_N + Mm_N
    minorA_N = 2*mm_N + Mm_N
    O_N = 2*OO_N

    percentCalled = '%.2f' % (100.0*(majorA_N + minorA_N)/(majorA_N + minorA_N + O_N)) if majorA_N + minorA_N + O_N > 0 else '0.00'

    if majorA == minorA  or sum(HW1(MM_N, Mm_N,mm_N)==0):
        p = 1.0
    else:
        hw, p = sc.chisquare(np.array([MM_N, Mm_N,mm_N]), f_exp=HW1(MM_N, Mm_N,mm_N), ddof=1)
        p = '%.5f' % p
    indel_N = N*4 - (majorA_N + minorA_N + O_N)

    if majorA != minorA:
        MAF = '%.4f' % (1.0*minorA_N / (majorA_N + minorA_N))
        return [majorA, minorA, MM_N, Mm_N, mm_N, OO_N, majorA_N, minorA_N, O_N, percentCalled,MAF,p]
    else:
        MAF = 0.0
    return [majorA, '', MM_N, Mm_N, mm_N, OO_N, majorA_N, minorA_N, O_N, percentCalled,MAF,p]

for k in range(statsP.shape[0]):
    CNT = verify(statsP[k], k)
    if  CNT !=  N*2:
        print >>sys.stderr, k, statsP[k], CNT, "!=", N*2, parents

hdb_stat = open('STAT/'+chunkN+'-stat.stat', 'wb')           
np.save(hdb_stat, np.vstack((statsP,statsC)))
hdb_stat.close()

print >>sys.stderr, "Done"

    
