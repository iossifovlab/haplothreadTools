#!/usr/bin/env python

import sys
from collections import Counter, defaultdict
import numpy as np
from scipy import stats

#if len(sys.argv) < 2:
#    print "Usage: expl.py <HDB>"
#    sys.exit(1) 

inFn = "D3.txt"

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

mFn = 'multiplex.txt'
qFn = 'quads.txt'
M = []
Q = []
D = [M,Q]
for i,f in enumerate([mFn,qFn]):
    with open(f) as F:
        for l in F:
            D[i].append(l.strip("\n\r"))


PF = open(inFn)
CNT = defaultdict()
for l in PF:
    cs = l.strip("\n\r").split("\t")
    fmId,prId,fId,mId,gender,aff = cs[:6]
    CNT[(fmId, prId, mId, gender, aff)] = int(cs[-1])
PF.close()

prb_m = [v for k,v in list(CNT.items()) if k[4] == '1' and k[3] == '1']
prb_f = [v for k,v in list(CNT.items()) if k[4] == '1' and k[3] == '2']
sib_m = [v for k,v in list(CNT.items()) if k[4] == '1' and k[3] == '1' and k[2] == '0']
sib_f = [v for k,v in list(CNT.items()) if k[4] == '1' and k[3] == '2' and k[2] == '0']

prb_M = [v for k,v in list(CNT.items()) if k[4] == '1' and k[3] == '1' and k[0] in M]
prb_Q = [v for k,v in list(CNT.items()) if k[4] == '1' and k[3] == '1' and k[0] in Q]
prb_Qf = [v for k,v in list(CNT.items()) if k[4] == '1' and k[3] == '2' and k[0] in Q]



