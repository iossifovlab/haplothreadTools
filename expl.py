#!/usr/bin/env python

import sys
from collections import Counter
import numpy as np

#if len(sys.argv) < 2:
#    print "Usage: expl.py <HDB>"
#    sys.exit(1) 

bimFn = "chip-sorted.bim"
pedFn = "SPARK-sorted.ped"

bim = []
with open(bimFn, 'r') as f:
    for l in f:
        cs = l.strip("\n\r").split("\t")
        bim.append([cs[0],cs[3],cs[4],cs[5]])

def processP(P):
    L = len(P[0])
    N = len(P)
    data = []
    for n in range(L-1):
        m = bim[n][2]
        M = bim[n][3]
        m1 = bim[n+1][2]
        M1 = bim[n+1][3]
        Z = '00'
        K = [M+M, M+m, m+M, m+m, Z]
        K1 = [M1+M1, M1+m1, m1+M1, m1+m1, Z]
        keys = [(x,y) for x in K for y in K1]
        cnt = Counter([(P[i][n],P[i][n+1]) for i in range(N)])
        data.append([cnt[k] for k in keys])    
    return np.array(data)

PF = open(pedFn)
P = []
n = 0
CNT = None
start = True
for l in PF:
    cs = l.strip("\n\r").split("\t")
    fmId,prId,fId,mId,gender,aff = cs[:6]
    if fId == '0' and mId == '0':
        P.append([cs[6+2*k]+cs[6+2*k+1] for k in range(len(cs[6:])/2)])
        n += 1
    if n % 1000 == 0 and len(P) > 0:
        if start:
            CNT = processP(P)
            start = False
        else:
            CNT += processP(P)
        P = []
        print >>sys.stderr, "n", n
    #if n == 2300:
    #    break

PF.close()

CNT += processP(P)

#"""
F = open("D.txt", 'w')
for d in CNT:
    F.write('\t'.join(map(str,d)) +"\n")
F.close()
#"""
