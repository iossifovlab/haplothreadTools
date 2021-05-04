#!/usr/bin/env python
import os,sys
#print 'hi'
from collections import defaultdict, Counter
from numpy import genfromtxt

DIR = 'SPARK_Freeze_Three/'
PSF = open(DIR+"array.genotype/sp_ids_in_genotype_ped.txt")
spIdsInPed = {l.strip("\n\r") for l in PSF}
PSF.close()

DT = genfromtxt(DIR + "SPARK.30K.Release.mastertable.20181105.txt", delimiter='\t',dtype=None,names=True, case_sensitive=True, encoding='ASCII')

print("There are", len({R['sf_id'] for R in DT}), "families", file=sys.stderr)

rlC = Counter([R['role'] for R in DT])

for rl,C in sorted(list(rlC.items()),key=lambda x: x[1]):
    print("%5d\t%s" % (C,rl), file=sys.stderr)

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


cleanFMS = {fId:RS for fId,RS in list(FMS.items()) if cleanF(RS)}
print("There are", len(cleanFMS), "clean families", file=sys.stderr)

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
for fmI,RS  in list(cleanFMS.items()):
    okFmSum[classifyOKFam(RS)][fmI] = RS

OF = open('cleanFamilyStats.txt','w')
hcs = "nAffectedParents nUnaffctedParents nAfftedChildren nUnaffctedChildren numberOfFamilies".split(" ")
#print "\t".join(hcs) 
OF.write("\t".join(hcs) + "\n")
for tp,c in sorted(list(okFmSum.items()),key=lambda x: len(x[1]),reverse=True):
    cs = list(map(str,list(tp) + [len(c)]))
    #print "\t".join(cs)
    OF.write("\t".join(cs) + "\n")
OF.close()

def main():
    quads = {f for f in cleanFMS if classifyOKFam(cleanFMS[f]) == (0,2,1,1)}

    with open('quads.txt', 'w') as f:
        for q in quads:
            f.write(q+'\n')


    multiplex = {f for f in cleanFMS if classifyOKFam(cleanFMS[f])[0] == 0 and classifyOKFam(cleanFMS[f])[1] == 2 and classifyOKFam(cleanFMS[f])[2] > 1}

    with open('multiplex.txt', 'w') as f:
        for q in multiplex:
            f.write(q+'\n')


    trios = {f for f in cleanFMS if classifyOKFam(cleanFMS[f])[0] == 0 and classifyOKFam(cleanFMS[f])[1] == 2 and classifyOKFam(cleanFMS[f])[2] == 1 and classifyOKFam(cleanFMS[f])[3] == 0}


    with open('trios.txt', 'w') as f:
        for q in trios:
            f.write(q+'\n')

if __name__ == "__main__":
	main()
