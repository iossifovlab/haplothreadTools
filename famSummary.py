#!/usr/bin/env python
import os,sys
#print 'hi'
from collections import defaultdict, Counter
from numpy import genfromtxt

DIR = os.getenv('PROJECT_DIR')
PSF = open(DIR+"/sp_ids_in_genotype_ped.txt")
spIdsInPed = {l.strip("\n\r") for l in PSF}
PSF.close()

DT = genfromtxt(DIR + "/mastertable.txt", delimiter='\t',dtype=None,names=True, case_sensitive=True)

print >>sys.stderr, "There are", len({R['sf_id'] for R in DT}), "families"

rlC = Counter([R['role'] for R in DT])

for rl,C in sorted(rlC.items(),key=lambda x: x[1]):
    print >>sys.stderr, "%5d\t%s" % (C,rl)

FMS = defaultdict(list)
for R in DT:
    FMS[R['sf_id']].append(R)

okRoles = set(['Father', 'Mother', 'Proband', 'Older Sibling', 'Younger Sibling'])

def cleanF(RS):
    nMoms = len([R for R in RS if R['role'] == 'Mother'])
    nDads = len([R for R in RS if R['role'] == 'Father'])
    nWeirdRls = len([R for R in RS if R['role'] not in okRoles])
    nUndiagnosed = len([R for R in RS if R['asd_flag']==0])
    nNoGenotype = len([R for R in RS if R['sp_id'] not in spIdsInPed])

    return nMoms == 1 and nDads == 1 and nWeirdRls == 0 and \
           len(RS) > 2 and nUndiagnosed == 0 and nNoGenotype == 0


cleanFMS = {fId:RS for fId,RS in FMS.items() if cleanF(RS)}
print >>sys.stderr, "There are", len(cleanFMS), "clean families"
import pickle

#sp_id	sf_id	father_id	mother_id	gender	asd_flag	role	age(year)	wes.CRAM	wes.GVCF	array.idat	array.signal_file
HD = {v:i for i,v in enumerate("sp_id sf_id father_id mother_id gender asd_flag role".split(" "))}
md = defaultdict()
for fId, R in cleanFMS.items():
    for p in R:
        s=1
        pId = p[HD['sp_id']]
        assert pId not in md
        dad = p[HD['father_id']]
        mom = p[HD['mother_id']]
        gender=p[HD['gender']]
        asd_flag=p[HD['asd_flag']]
        role=p[HD['role']]
        md[pId] = {'fam':fId, 'indId':pId, 'dad': dad, 'mom':mom, 'gender':gender,
                   'role':role, 'aff':asd_flag}

with open('fms.pckl', 'wb') as f:
    pickle.dump(md, f)

def classifyOKFam(RS):
    pars = [R for R in RS if R['role'] in ['Mother','Father']]
    prbs = [R for R in RS if R['role'] == 'Proband'] 
    sibs = [R for R in RS if R['role'] not in ['Mother','Father', 'Proband'] ]
    
    def au(RS):
        csts = Counter([R['asd_flag'] for R in RS])
        assert set(csts.keys()).issubset(set([1, 2]))
        return [csts[2], csts[1]]
    # return tuple(au(pars) + au(prbs) + au(sibs) + au(prbs+sibs))
    return tuple(au(pars) + au(prbs+sibs))

okFmSum = defaultdict(dict)
for fmI,RS  in cleanFMS.items():
    okFmSum[classifyOKFam(RS)][fmI] = RS

OF = open('cleanFamilyStats.txt','w')
hcs = "nAffectedParents nUnaffctedParents nAfftedChildren nUnaffctedChildren numberOfFamilies".split(" ")
#print "\t".join(hcs) 
OF.write("\t".join(hcs) + "\n")
for tp,c in sorted(okFmSum.items(),key=lambda x: len(x[1]),reverse=True):
    cs = map(str,list(tp) + [len(c)])
    #print "\t".join(cs)
    OF.write("\t".join(cs) + "\n")
OF.close()


quads = {f for f in cleanFMS if classifyOKFam(cleanFMS[f]) == (0,2,1,1)}
#"""
with open('quads.txt', 'w') as f:
    for q in sorted(quads):
        f.write(q+'\n')
#"""

multiplex = {f for f in cleanFMS if classifyOKFam(cleanFMS[f])[0] == 0 and classifyOKFam(cleanFMS[f])[1] == 2 and classifyOKFam(cleanFMS[f])[2] > 1}
#"""
with open('multiplex.txt', 'w') as f:
    for q in sorted(multiplex):
        f.write(q+'\n')
#"""

trios = {f for f in cleanFMS if classifyOKFam(cleanFMS[f])[0] == 0 and classifyOKFam(cleanFMS[f])[1] == 2 and classifyOKFam(cleanFMS[f])[2] == 1 and classifyOKFam(cleanFMS[f])[3] == 0}

#"""
with open('trios.txt', 'w') as f:
    for q in sorted(trios):
        f.write(q+'\n')

#"""
