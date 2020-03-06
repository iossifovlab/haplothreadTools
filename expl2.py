#!/usr/bin/env python

import sys
from collections import Counter
import numpy as np

#if len(sys.argv) < 2:
#    print "Usage: expl.py <HDB>"
#    sys.exit(1) 

bimFn = "chip-sorted.bim"
pedFn = "SPARK-sorted.ped"

"""
bim = []
with open(bimFn, 'r') as f:
    for l in f:
        cs = l.strip("\n\r").split("\t")
        bim.append([cs[0],cs[3],cs[4],cs[5]])
"""

keys = ['00',
 'AA',
 'AC',
 'AG',
 'AT',
 'CA',
 'CC',
 'CG',
 'CT',
 'GA',
 'GC',
 'GG',
 'GT',
 'TA',
 'TC',
 'TG',
 'TT']

hom = ['AA','CC','GG','TT']

PF = open(pedFn)
F = open("D2.txt", 'w')
n = 0
CNT = []
start = True
for l in PF:
    cs = l.strip("\n\r").split("\t")
    fmId,prId,fId,mId,gender,aff = cs[:6]
    D = Counter([cs[6+2*k]+cs[6+2*k+1] for k in range(len(cs[6:])/2)])
    H = sum([D[k] for k in hom])
    F.write('\t'.join(map(str,[D[k] for k in keys] + [H])) +"\n")
    if n % 1000 == 0:
        print >>sys.stderr, "n", n
    n += 1
    #if n == 10:
    #    break

PF.close()
F.close()



