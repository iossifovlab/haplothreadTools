#!/usr/bin/env python
import os, sys
from collections import defaultdict, Counter
from glob import glob
import numpy as np

fn = sorted(glob('STATS/posStats-*.npz'))
print(fn, file=sys.stderr)

with open(fn[0], 'rb') as F:
        data = np.load(F, encoding='ASCII')
        data.allow_pickle=True
        statsPQ = data.get('SPQ')
        statsPT = data.get('SPT')
        statsPM = data.get('SPM')        
        statsCQ = data.get('SCQ')
        statsCT = data.get('SCT')
        statsCM = data.get('SCM')        
        parentsQ = data.get('pq')
        parentsT = data.get('pt')
        parentsM = data.get('pm')        
        print ('loaded', fn[0], file=sys.stderr)
        
for f in fn[1:]:
    with open(f, 'rb') as F:
        data = np.load(F, encoding='ASCII')
        data.allow_pickle=True
        statsPQ += data.get('SPQ')
        statsPT += data.get('SPT')
        statsPM += data.get('SPM')        
        statsCQ += data.get('SCQ')
        statsCT += data.get('SCT')
        statsCM += data.get('SCM')        
        parentsQ += data.get('pq')
        parentsT += data.get('pt')
        parentsM += data.get('pm')        
        print ('loaded', f, file=sys.stderr)
        
with open('statsSummaryQ.npz', 'wb') as F:
    np.savez(F, SPQ=statsPQ, SCQ=statsCQ, pq=parentsQ)

with open('statsSummaryT.npz', 'wb') as F:
    np.savez(F, SPT=statsPT, SCT=statsCT, pt=parentsT)

with open('statsSummaryM.npz', 'wb') as F:
    np.savez(F, SPM=statsPM, SCM=statsCM, pm=parentsM)

