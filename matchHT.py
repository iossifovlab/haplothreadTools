#!/bin/env python
from math import *
from pylab import *
import pylab as plt
import numpy as np
import scipy.stats as sci
import glob 
from time import time
from collections import Counter
interactive(1)
import sys
from time import time

if len(sys.argv) < 2:
    print("Usage: statsHDB.py <HDB>")
    sys.exit(1) 

HDB = sys.argv[1]

expectedCodes = "ACGTMRWSYKcdou"
MAP = {v:i for i,v in enumerate(expectedCodes)}
RMAP = {v:k for k,v in list(MAP.items())}
ambiguous = {
    'AC':'M', 'CA':'M', 
    'AG':'R', 'GA':'R', 
    'AT':'W', 'TA':'W', 
    'CG':'S', 'GC':'S', 
    'CT':'Y', 'TC':'Y',
    'GT':'K', 'TG':'K'
}

def loadHDB(HDB):
    HT = {}
    HF = open(HDB + "-hpth.txt")
    HF.readline()
    for l in HF:
        fmId,prId,hpId,hpth = l.strip("\n\r").split("\t")
        HT[(fmId,prId,hpId)] = np.array([MAP[i] for i in hpth], dtype=int8)
    HF.close()

    POS = []
    n = -1
    HF = open(HDB + "-pos.txt")
    HF.readline()
    for l in HF:
        n += 1
        POS.append(l.strip("\n\r").split("\t"))
    HF.close()

    pos = np.array([k[1] for k in POS]).astype(int32)
    return [HT,pos]

# mask is 0 for possible match and 1 for mismatch

mask = {(x,y):0 for x in expectedCodes for y in expectedCodes}
for x in 'ACGTMRWSYK':
    for y in 'ACGTMRWSYK':
        mask[(x,y)] = 1
        mask[(y,x)] = 1

for x in "ACGTMRWSYK":
    mask[(x,x)] = 0

for k,v in list(ambiguous.items()):
    mask[(k[0],v)] = 0
    mask[(k[1],v)] = 0
    mask[(v,k[0])] = 0
    mask[(v,k[1])] = 0

for P in ['MRW', 'RSK', 'MSY', 'WYK']:
    for a in P:
        for b in P:
            mask[(a,b)] = 0
            mask[(b,a)] = 0

M = np.array([mask[(x,y)] for x in expectedCodes for y in expectedCodes]).reshape(len(expectedCodes), len(expectedCodes))

def compareHT(HT, pos, k1,k2):
    """ takes dictionary of haplothreads (HT), positions (pos) and two keys (k1,k2) into HT
    and returns list of mathshing intervals marked if they have detected matches inside or not
    and a list of length also marked
    """

    C = M[HT[k1], HT[k2]]
    id =np.where(np.array(C) == 1)[0]

    ints = [[id[i-1],id[i],id[i]-id[i-1]-1] for i in range(1,len(id))]
    if len(id) == 0:
        return [[0,len(C)-1], [pos[0], pos[-1]], [pos[-1] - pos[0]], [pos[-1]-pos[0]]]
    if id[0] > 0:
        ints.insert(0,[0,id[0], id[0]-1])
    if id[-1] != len(C)-1:
        ints.append([id[-1],len(C)-1, len(C)-id[-1]-2])
    ints = np.array(ints)
    gints = np.array([[pos[k[0]], pos[k[1]], k[2]] for k in ints])
    L = gints[:,1] - gints[:,0]
    id0 = np.where(gints[:,2] > 0)
    L0 = gints[id0,1] - gints[id0,0]
    return [ints, gints, L, L0[0]]

N = 680658

HT, pos = loadHDB(HDB)
for n in sorted(np.random.randint(1,N, 100)):
    print('\t'.join(map(str,[22, pos[n], pos[n]])))





"""
HT, pos = loadHDB(HDB)
keys = HT.keys()

ints, gints, L, L0 = compareHT2(HT, pos, keys[3], keys[7])


for k1 in range(7):
    for k2 in range(k1+1,8):
        C = M[HT[keys[k1]], HT[keys[k2]]]
        D = C[1:] - C[:-1]
        a = [i for i in range(len(D)-1) if D[i] == 1 and D[i+1] == 1]
        print len(a)
"""
