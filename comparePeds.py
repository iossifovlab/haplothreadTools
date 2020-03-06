#!/usr/bin/env python
import sys,os
import numpy as np
from itertools import izip
from collections import defaultdict, Counter
import pickle

if len(sys.argv) < 6:
    print "Usage: comparePeds.py <ped1> <ped2> <bim1> <bim2> <posFn>"
    exit()

ped1Fn = sys.argv[1]
ped2Fn = sys.argv[2]
bim1Fn = sys.argv[3]
bim2Fn = sys.argv[4]
posFn = sys.argv[5]

refA = []
with open(posFn, 'r') as f:
    f.readline()
    for l in f:
        cs = l.strip("\n\r").split("\t")
        refA.append(cs[2])

with open(bim1Fn, 'r') as f1, open(bim2Fn, 'r') as f2:
    i = 0
    for l1,l2 in izip(f1,f2):
        assert (l1 and l2 or (not l1 and not l2))
        if not l1: break
        cs1 = l1.strip("\n\r").split("\t")
        cs2 = l2.strip("\n\r").split("\t")        
        assert (cs1[0],cs1[3]) == (cs2[0],cs2[3])
        try:
            if cs1[4] == '0':
                assert cs1[5] in set([cs2[4],cs2[5]])
            else:
                assert set([cs1[4],cs1[5]]) == set([cs2[4],cs2[5]])
            assert cs2[4] == refA[i] or cs2[5] == refA[i]
        except:
            print >>sys.stderr, "bim alleles do not match", cs1, cs2
        i += 1

def T(x):
    #return 'het' if x[1] != x[2] else 'homR' if x[0] == x[1] else '0' if x[1] == '0' else 'homA'
    return 1 if x[1] != x[2] else 0 if x[0] == x[1] else 3 if x[1] == '0' else 2

TT={}
for i in range(4):
    for k in range(4):
        TT[(i,k)]=[0 for n in range(16)]
        TT[(i,k)][4*i+k] = 1

    
L = len(refA)
n = 0

start = True
with open(ped1Fn, 'r') as f1, open(ped2Fn, 'r') as f2:
    for l1,l2 in izip(f1,f2):
        assert (l1 and l2 or (not l1 and not l2))
        if not l1: break

        cs1 = l1.strip("\n\r").split("\t")
        cs2 = l2.strip("\n\r").split("\t")

        assert len(cs1) == len(cs2)
        assert len(cs1[6:]) == L*2
        assert (tuple(cs1[:6]) == tuple(cs2[:6]))

        g1 = [T([refA[k], cs1[6+2*k], cs1[7+2*k]]) for  k in range(L)]
              
        g2 = [T([refA[k], cs2[6+2*k], cs2[7+2*k]]) for  k in range(L)]
        cnt = np.sum(np.array([TT[(x,y)] for x,y in zip(g1,g2)]), axis=0)
        if start:
            CNT = cnt
            start= False
        else:
            CNT += cnt
        n += 1
        #print n
        if n % 100 == 0:
            print >>sys.stderr, "processed", n, "lines"
        #    break

h_cnt = np.zeros((4,4))
hMap={x:y for x,y in zip('homR het homA 0'.split(' '), range(4))}
hMapr = {y:x for x,y in hMap.items()}

h_cnt=CNT.reshape((4,4))

S1 = np.sum(h_cnt, axis = 0)
S2 = np.sum(h_cnt, axis = 1)
S = sum(S1)
assert (S == sum(S2))

print '\t'.join('* homR het homA 0 Sum'.split(' '))
for k in range(4):
    out = [hMapr[k]] + list(h_cnt[k,:]) + [S2[k]]
    print '\t'.join(map(str,out))
print '\t'.join(map(str,['Sum'] + list(S1) + [S]))
