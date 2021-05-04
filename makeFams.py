#!/usr/bin/env python
import numpy as np
from time import time
import sys,os
from collections import defaultdict

fn = os.environ['SPARKMASTERTABLE']

pFn = 'persons.txt'

persons = {}
n = 0
with open(pFn, 'r') as f:
    for l in f:
        persons[l.strip('\n\r')] = n
        n += 1

hd = 'sp_id sf_id father_id mother_id gender asd_flag role age(year)'.split(' ')

fam = defaultdict(list)
quads = defaultdict(list)

class P:
    pass

with open(fn, 'r') as f:
    HD = {h:i for i,h in enumerate(f.readline().strip('\n\r').split('\t'))}
    for l in f:
        cs = l.strip('\n\r').split('\t')
        rec = {h:cs[HD[h]] for h in hd}
        p = P()
        p.famId = rec['sf_id']
        p.id = rec['sp_id']
        p.dad = rec['father_id']
        p.mom = rec['mother_id']
        p.role = rec['role']
        p.gender = rec['gender']
        p.age = rec['age(year)']
        p.aff = rec['asd_flag']
        p.n = persons[p.id] if p.id in persons else -1
        fam[rec['sf_id']].append(p)

def roleN(f, role):
    return sum([1 for p in f if p.role == role])

for fId,f in list(fam.items()):
    momN = roleN(f,'Mother')
    dadN = roleN(f, 'Father')
    unAffParents = sum([1 for p in f if (p.role == 'Father' and p.aff == '1')
                        or (p.role == 'Mother' and p.aff == '1')])
    prbN = sum([1 for p in f if p.role == 'Proband' and p.aff == '2'])

    unAffChildren = sum([1 for p in f if (p.aff == '1' and 'Sibling' in p.role)])
    
    weirdN = sum([1 for p in f if not (p.role == 'Mother'
                                       or p.role == 'Father'
                                       or p.role == 'Proband'
                                       or p.role == 'Older Sibling'
                                       or p.role == 'Younger Sibling')])
    N = len(f)
    genN = sum([1 for p in f if p.n > -1])
    if (momN == 1 and
        dadN == 1 and
        weirdN == 0 and
        unAffParents == 2 and
        unAffChildren == 1 and
        N == 4 and
        prbN == 1 and
        genN == 4):
        quads[fId] = f
"""
for q in quads:
    print(q)
"""
