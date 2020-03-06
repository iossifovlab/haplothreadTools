#!/usr/bin/env python
from math import *
import numpy as np
import sys
from collections import defaultdict, Counter
from time import time
# import matplotlib as mpl
# mpl.use('Agg')

L = 10000000
# L = 1000
if len(sys.argv) <3:
    print "Usage: compare.py <hdb1> <hdb2>"
    sys.exit(1)


#HT1,HT2 = "NYGC-chip551839-SSCquads-good372556", "NYGC-chip551839-SSCquads-good372556-WGS"

HT1,HT2 = sys.argv[1:3]

def load_htdb(hdbName):
    #print("loading", hdbName, "...",file=sys.stderr)
    PSs = []
    HTs = []
    with open(hdbName + '-pos.txt', 'r') as F:
        F.readline()
        for l in F:
            ch,ps,rf,A,C,G,T = l.strip('\n\r').split('\t')
            ps,A,C,G,T = map(int,[ps,A,C,G,T])
            PSs.append((ch,ps,rf,A,C,G,T))
    with open(hdbName + '-hpth.txt', 'r') as F:
        F.readline()
        for l in F:
            fId,pId,hId,hth = l.strip('\n\r').split('\t')
            HTs.append((fId,pId,hId,hth[:L]))
    return PSs,HTs

PSs1,HTs1 = load_htdb(HT1)
PSs2,HTs2 = load_htdb(HT2)

# for c in [0,1,2]:
for c in [0,1]:
    assert [x[c] for x in PSs1] == [x[c] for x in PSs2]

for c in [0,1,2]:
    assert [x[c] for x in HTs1] == [x[c] for x in HTs2]

unknownCs = set('dcuo')
rvcmpD = [
    'ACGTRYSWKM',
    'TGCAYRSWMK'
]
rvcmpM = {}
for f,t in zip(rvcmpD[0],rvcmpD[1]):
    rvcmpM[f] = t
for u in unknownCs:
    rvcmpM[u] = u
def revCmpV(v):
    return [rvcmpM[n] for n in v]
    
TT = []
# for p in range(0,len(PSs1),100):
nBad = 0
for p in range(0,len(HTs1[0][3])):
    c1 = [h[3][p] for h in HTs1]
    c2 = [h[3][p] for h in HTs2]

    def smry(a,b):
        if a in unknownCs and b in unknownCs:
            return 'b'
        if a in unknownCs:
            return 'f'
        if b in unknownCs:
            return 's'
        if a == b:
            return 'M'
        else:
            return 'n'


    def match(c1,c2):
        cnts = Counter([smry(a,b) for a,b in zip(c1,c2)])
        M = cnts['M']
        n = cnts['n']
        if n+M:
            n_p = float(n) / (n+M)
        else:
            n_p = -1.0 
        return n_p,n,M
    flp = 'N'
    n_p,n,M = match(c1,c2)

    if n_p > 0.5:
        flp = 'Y'
        n_p,n,M = match(revCmpV(c1),c2)

    TT.append((PSs1[p][:2],n_p,M+n))
    print "\t".join(map(str,[PSs1[p][0],PSs1[p][1],n_p,n,M+n,flp]))
    # if n_p > 0.01:
    # if n_p > 0.003 and n_p < 0.01:
    #     nBad +=1 
    #     cntsR = Counter([(a,b) for a,b in zip(c1,c2)])
    #     print(nBad,PSs1[p][:2],n_p,M+n,cnts)
    #     print("\t",cntsR)
    #     print("\t",PSs1[p])
    #     print("\t",PSs2[p])

'''    
clf()
subplot(3,1,1)
plot([x[2] for x in TT])
a = gca()
subplot(3,1,2,sharex=a)
plot([x[1] for x in TT],'.')
subplot(3,1,3,sharex=a)
plot([x[0][1] for x in TT])
show()
# savefig('ivan.png')
'''
