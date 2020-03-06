#!/bin/env python
from math import *
import numpy as np
import sys

if len(sys.argv) < 4:
    print "Usage: statsHDB.py <HDB> <numPos> <numHT>"
    sys.exit(1) 

HDB = sys.argv[1]
numPos = int(sys.argv[2])
numHT = int(sys.argv[3])
N = 100
X = np.random.choice(range(numPos), N, replace=False)
Y = np.random.choice(range(numHT), N, replace=False)
assert (len(set(X)) == N)
assert (len(set(Y)) == N)
id = np.argsort(Y)
Y = Y[id]
X = X[id]

#print X
#print Y
#print HDB
HT = {}
def createRP(HDB, X, Y):
    HF = open(HDB + "-hpth.txt")
    HF.readline()
    n = -1
    k = 0
    for l in HF:
        n += 1
        if n != Y[k]:
            continue
        fmId,prId,hpId,hpth = l.strip("\n\r").split("\t")
        HT[(fmId,prId,hpId,X[k])] = [hpth[X[k]]]
        k += 1
        #print k
        if k == N:
            break
    HF.close()

    pos2ht = { k[3]:k for k in HT.keys()}
    X = sorted(X)
    HF =  open(HDB + "-pos.txt")
    HF.readline()
    n = -1
    k = 0
    for l in HF:
        n += 1
        if n != X[k]:
            continue
        cs = l.strip("\n\r").split("\t")
        HT[pos2ht[X[k]]].append(cs[1])
        k += 1
        if k == N:
            break
    HF.close()

createRP(HDB, X, Y)


for k in sorted(HT):
    print '\t'.join(map(str, list(k) + HT[k]))





