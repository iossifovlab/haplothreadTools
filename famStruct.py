#!/usr/bin/env python

import sys, os, gzip
import numpy as np
from famSummary import *


#quads = 'goodFam-noMZ.txt'
#quadsFn = 'quads.txt'
#quadsFn = 'multiplex.txt'

if len(sys.argv) <2:
    print "Usage: famStruct.py <file with the list of families>"
    exit()

famType = sys.argv[1]
if famType == 'Q':
    fams = quads
elif famType == 'T':
    fams = trios
elif famType == 'M':
    fams = multiplex
else:
    print "wrong famType", famType, "should be on of Q, T, M"
    exit()

#with open(quadsFn, 'r') as f:
#    for l in f:
#        quadFams.append(l.strip('\n\r'))

GEN = {'1':'M','2':'F', 1:'M', 2:'F'}
AFF = {'1':'sib', '2':'prb', 1:'sib',2:'prb'}
ROLE = {'Mother':'mom', 'Father':'dad', 'Older Sibling':'sib1', 'Younger Sibling':'sib2', 'Proband':'prb'}
persons = {}

famCnt = {f:defaultdict(int) for f in fams}
for f in fams:
    fam = cleanFMS[f]
    for p in fam:
        if ROLE[p[6]] == 'mom':
            pId = p[1]+'.mo'
            famCnt[f]['mom'] += 1
        elif ROLE[p[6]] == 'dad':
            pId = p[1]+'.fa'
            famCnt[f]['dad'] += 1
        elif AFF[p[5]] == 'prb':
            famCnt[f]['prb'] += 1
            pId = p[1]+'.p' + str(famCnt[f]['prb'])
        elif AFF[p[5]] == 'sib':
            famCnt[f]['sib'] += 1            
            pId = p[1]+'.s' + str(famCnt[f]['sib'])
        persons[p[0]] = [p[0], p[1], pId,GEN[p[4]]]

for p,v in sorted(persons.items(), key=lambda x: [x[1][1],x[1][2]]):
    print '\t'.join(persons[p])


