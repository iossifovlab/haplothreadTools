#!/usr/bin/env python
import os, sys
#from famSummary import *
from collections import defaultdict, Counter
from glob import glob
from hw import HW1
import numpy as np
import scipy.stats as sc
from famSummary import *


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

if len(sys.argv) <2:
    print("Usafe statsSummary2.py <letter Q|T|M>")
    exit()

famType = sys.argv[1]
MAP = []
with open('SPARK.30K.genotype.release_set.map', 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        MAP.append([cs[0], cs[3]])


#N = 5802*2

with open('statsSummary' + famType + '.npz', 'rb') as F:
    #statsP, statsC, parents = np.load(F, encoding='ASCII')
    data = np.load(F, encoding='ASCII')
    data.allow_pickle = True
    statsP  = data.get('SP'+famType)
    statsC  = data.get('SC'+famType)
    parents = data.get('p'+ famType.lower())
    
N = parents
#"""
def verify(sP, k):
    basesP = [item for sublist in sP for item in sublist]
    c = {x:sum([sP[y]*2 if (x,x) == y else sP[y]
                for y in sP if x in y]) for x in basesP}
    return sum(c.values())

M={v:i for i,v in enumerate('0TGCA')}
def m(x):
    return M[x] if x in M else 5 
    
def alleles(sP, sC, k):
    basesP = [item for sublist in sP for item in sublist]
    c = {x:sum([sP[y]*2 if (x,x) == y else sP[y]
                for y in sP if x in y]) for x in basesP}
    d = sorted(list(c.items()), key=lambda x: (x[1],m(x[0])), reverse=True)
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
        d = sorted(list(c.items()), key=lambda x: x[1], reverse=True)
        print('\t'.join(map(str,[MAP[k][0],MAP[k][1],
                                         'parents', sP, 'children', sC, d])), file=sys.stderr)
        d = [x for x in d if x[0] != '0']
        
        if len(d) > 1:
            minorA = d[1][0]
        if len(d) > 2:
            print("Strange", file=sys.stderr)
    
    keys = [[(majorA, majorA)],
            [(minorA, minorA)],
            [(majorA, minorA), (minorA, majorA)],
            [('0','0')]]

    MM_N, mm_N, Mm_N, OO_N = [sum([sP[i] for i in k]) for k in keys]

    if minorA == majorA:
        mm_N = 0
        Mm_N = 0
    majorA_N = 2*MM_N + Mm_N
    minorA_N = 2*mm_N + Mm_N
    O_N = 2*OO_N

    percentCalled = '%.2f' % (100.0*(majorA_N + minorA_N)/(majorA_N + minorA_N + O_N))
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


#"""
for k in range(statsP.shape[0]):
    CNT = verify(statsP[k], k)
    if  CNT !=  N*2:
        print(k, statsP[k], CNT, "!=", N*2, parents, file=sys.stderr)

#"""
#"""        
out = []

print('\t'.join('chrom pos majorA minorA MM_N Mm_N mm_N 00_N majorA_N minorA_N 0_N percentCalled MAF HW_p'.split(' ')))
for k in range(statsP.shape[0]):
    out = MAP[k]
    out = out+alleles(statsP[k],statsC[k], k)
    print('\t'.join(map(str,out)))

               
#"""
