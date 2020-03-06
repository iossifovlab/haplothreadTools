#!/usr/bin/env python

import sys
from collections import Counter
import numpy as np

if len(sys.argv) < 2:
    print("Usage: statsHDB.py <HDB>")
    sys.exit(1) 

HDB = sys.argv[1]

posA = []
PF = open(HDB + "-pos.txt")
PF.readline()
for l in PF:
    class POS:
        pass

    p = POS()
    p.ch,p.ps,p.refA,p.A,p.C,p.G,p.T = l.strip("\n\r").split("\t")
    p.ps,p.A,p.C,p.G,p.T = list(map(int,[p.ps,p.A,p.C,p.G,p.T]))
    
    posA.append(p)
    
PF.close()

refHTA = np.array([p.refA for p in posA])

expectedCodes = "ACGTcdouMRWSYK"
expectedCodesSet = set(expectedCodes)

print("\t".join("familyId personId haplotypeId".split() + list(expectedCodes) + ["nRef"]))
HF = open(HDB + "-hpth.txt")
HF.readline()
for l in HF:
    fmId,prId,hpId,hpth = l.strip("\n\r").split("\t")
    cnts = Counter(hpth)
    assert not set(cnts) - expectedCodesSet
   
    htA = np.array(list(hpth)) 
    nRef = (htA==refHTA).sum()
    print("\t".join([fmId,prId,hpId] + \
            [str(cnts[c]) for c in expectedCodes] + [str(nRef)]))
HF.close()
