#!/usr/bin/env python
import sys,os
from collections import defaultdict
import pickle

if len(sys.argv) < 3:
    print("Usage: createPedMap.py <ped file> <pedigree file>")
    exit(1)

pedFn=sys.argv[1]
#fmsFn=sys.argv[2]
pedigreeFn = sys.argv[2]
#md=defaultdict()
H='fam indId dad mom gender role'.split(' ')

"""
F = open(fmsFn)
families,badFamilies = pickle.load(F)
F.close()
"""

#NYGCSampleID	sfid	spid	father	mother	sex	asd	role
md = defaultdict()
with open(pedigreeFn) as f:
    hd = {v:i for i,v in enumerate(f.readline().strip("\n\r").split("\t"))}
    for l in f:
        cs = l.strip("\n\r").split("\t")
        if "B01" in os.environ["PROJECT_DIR"] or "B02" in os.environ["PROJECT_DIR"]:
            md[cs[hd["NYGCSampleID"]]] = cs
        else:
            md[cs[hd["spid"]]] = cs

"""    
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
"""
n = 0
with open(pedFn, 'r') as f:
    # #N sampleId
    for l in f:
        n += 1
        cs = l.strip('\n\r').split('\t')
        if cs[1] in md:
            out = []
            rec = md[cs[1]]
            out.append(rec[hd["sfid"]])
            out.append(rec[hd["spid"]])
            out.append(rec[hd["father"]])
            out.append(rec[hd["mother"]])
            out.append(rec[hd["sex"]])
            out.append(rec[hd["asd"]])
            print('\t'.join(out + cs[6:]))
        #if n == 30:
        #    break

 
 
