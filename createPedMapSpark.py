#!/usr/bin/env python
import sys,os
from collections import defaultdict, Counter
import pickle

if len(sys.argv) < 2:
    print("Usage: createPedMap.py <ped file>")
    exit()

pedFn=sys.argv[1]

F = open('fms.pckl')
md = pickle.load(F)
F.close()

"""
for fd in families.values():
    for pd in fd.memberInOrder:
        smId = os.path.basename(pd.atts['bamF']).split(".")[0]
        assert smId not in md
        famId = fd.familyId
        indId = pd.personId
        dad = famId + '.fa' if pd.role in ['sib', 'prb'] else '0'
        mom = famId + '.mo' if pd.role in ['sib', 'prb'] else '0'
        gender = pd.gender
        role = pd.role
        md[smId] = {'fam': famId, 'indId':indId, 'dad':dad, 
                    'mom':mom, 'gender':gender, 'role':role}
"""
gMap = {'F':'2', 'M':'1'}
aMap = {'prb':'2', 'sib':'1'}
H='fam indId dad mom gender aff'.split(' ')
with open(pedFn, 'r') as f:
    # #N sampleId
    for l in f:
        cs = l.strip('\n\r').split(' ')
        if cs[1] in md:
            out=md[cs[1]]
            print('\t'.join(list(map(str,[out[k] for k in H ]))+cs[6:]))
    

 
 
