#!/usr/bin/env python
import sys,os
from collections import defaultdict
import pickle

if len(sys.argv) < 3:
    print("Usage: createPedMap.py <ped file> <fms.pckl file for SSC2380>")
    exit(1)

pedFn=sys.argv[1]
fmsFn=sys.argv[2]
md=defaultdict()
H='fam indId dad mom gender role'.split(' ')

F = open(fmsFn)
families,badFamilies = pickle.load(F)
F.close()

for fd in list(families.values()):
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

gMap = {'F':'2', 'M':'1'}
aMap = {'prb':'2', 'sib':'1'}

with open(pedFn, 'r') as f:
    # #N sampleId
    for l in f:
        cs = l.strip('\n\r').split('\t')
        if cs[1] in md:
            out=md[cs[1]]
            famId = out['fam']
            out['gender']=gMap[out['gender']]
            out['affectStatus'] = '1' if out['role'] in ['mom','dad'] else aMap[out['role']]
            print('\t'.join([out[k] for k in H[:-1] 
                             + ['affectStatus']]+cs[6:]))
    

 
 
