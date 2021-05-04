#!/usr/bin/env python
import os, sys
from famSummary import *
from collections import defaultdict, Counter
import numpy as np

fams = list(cleanFMS.keys())

if len(sys.argv) < 3:
    print("Usage cCnt.py <chunk #> <ped file>")
    exit()

#pedFn = 'SPARK.30K.snp.SP.ped' 
chunk = sys.argv[1]
pedFn = sys.argv[2]
sline = int(chunk)*500
eline = (int(chunk)+1)*500


print(sline, eline, pedFn, "number of fams:", len(fams), file=sys.stderr)

ped = []
n=-1
counts = []
with open(pedFn, 'r') as f:
    for l in f:
        n +=1
        if n < sline:
            continue
        if n >= eline:
            break
        cs = l.strip('\n\r').split(' ')
        if not cs[0] in fams:
            continue
        pId = cs[1]
        G = cs[6:]
        counts.append([cs[0],cs[1],Counter(G)])
        #if n > 30:
        #    break

fn = open('STATS/cCnt-'+chunk+'.npz', 'wb')
np.savez(fn,counts=counts)
fn.close()

print("Done", 'STATS/cCnt-'+str(sline)+'-'+str(eline)+'.npz', file=sys.stderr) 
