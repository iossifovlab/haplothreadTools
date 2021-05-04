#!/usr/bin/env python
import numpy as np
from collections import defaultdict
import sys, os
from makeFams import *

if len(sys.argv) < 6:
    print("Usage: ped2HT.py <famId> <outDir> <filteredPos file> <ped file> <fitered families file>")
    exit()
    
famId = sys.argv[1]
outDir = sys.argv[2]
mapFn = sys.argv[3]
pedFn = sys.argv[4]
filteredFamsFn = sys.argv[5]

# fam comes from makeFams importan attribute is the person line number in ped
family = fam[famId]

filteredFams = defaultdict(list)
with open(filteredFamsFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        filteredFams[cs[0]].append(cs[1])
        
# read filteredPos file and store chrom, pos, refA and altA
#chrom   pos     majorA  minorA  MM_N    Mm_N    mm_N    00_N    majorA_N        minorA_N        0_N     percentCalled   MAF  HW_p
snpMap = []
with open(mapFn, 'r') as f:
    cs = f.readline().strip('\n\r').split('\t')
    # if no header use the first line
    if cs[0] == '1':
        snpMap.append(cs[:4])
    for l in f:
        cs = l.strip('\n\r').split('\t')
        snpMap.append(cs[:4])

family = [p for p in family if p.id in filteredFams[famId]]

sf = sorted([[p,int(p.n)] for p in family], key=lambda x: x[1])

nums = [s[1] for s in sf]
pers = [s[0] for s in sf]

FAM = {}

ROLE = {'Mother':'mom',
        'Father':'dad',
        'Proband':'prb',
        'Younger Sibling':'sib',
        'Older Sibling':'sib'}

# read family
n = 0
k = 0
children = []
with open(pedFn, 'r') as f:
    for l in f:
        # go throus all lines in ped file, but parse and process only onces
        # correspond to family member line numbers nums
        if n in nums:
            id = nums.index(n)
            role = ROLE[pers[id].role]
            cs = l.strip('\n\r').split(' ')
            fId, pId, faId, moId, sex, aff = cs[:6]
            if role == 'mom':
                mom = pId
            elif role == 'dad':
                dad = pId
            else:
                children.append(pId)
            G = cs[6:]
            assert (fId == famId)
            assert (pers[id].id == pId)
            try:
                assert (pers[id].gender == sex)
            except AssertionError:
                print("famId", fId, "pId", pId, 'sex mismatch: ', pers[id].gender, "!=", sex, file=sys.stderr)
            FAM[pId] = {'personId':pId,'genotype':np.array(G), 'role':role}
            k += 1
        n += 1
        if k == len(family):
            break

# end read family
members = [mom, dad] + sorted(children)

GS = np.vstack(tuple(FAM[p]['genotype'] for p in members))

ambiguous = {
    'AC':'M', 'CA':'M', 
    'AG':'R', 'GA':'R', 
    'AT':'W', 'TA':'W', 
    'CG':'S', 'GC':'S', 
    'CT':'Y', 'TC':'Y',
    'GT':'K', 'TG':'K'
}

gn2ph = {
    'cc cc cc':['altBase','altBase','altBase','altBase'],
    'cc ac cc':['altBase','altBase','altBase','refBase'],
    'cc ac ac':['altBase','altBase','refBase','altBase'],
    'cc aa ac':['altBase','altBase','refBase','refBase'],
    'ac cc cc':['altBase','refBase','altBase','altBase'],
    'ac cc ac':['refBase','altBase','altBase','altBase'],
    'ac ac cc':['altBase','refBase','altBase','refBase'],
    'ac ac ac':['E', 'E', 'E', 'E'],
    'ac ac aa':['refBase','altBase','refBase','altBase'],
    'ac aa ac':['altBase','refBase','refBase','refBase'],
    'ac aa aa':['refBase','altBase','refBase','refBase'],
    'aa cc ac':['refBase','refBase','altBase','altBase'],
    'aa ac ac':['refBase','refBase','altBase','refBase'],
    'aa ac aa':['refBase','refBase','refBase','altBase'],
    'aa aa aa':['refBase','refBase','refBase','refBase']
}

# handle X we have to augment this dict

gn2ph = {
    'cc cc cc':['altBase','altBase','altBase','altBase'],
    'cc ac cc':['altBase','altBase','altBase','refBase'],
    'cc ac ac':['altBase','altBase','refBase','altBase'],
    'cc aa ac':['altBase','altBase','refBase','refBase'],
    'ac cc cc':['altBase','refBase','altBase','altBase'],
    'ac cc ac':['refBase','altBase','altBase','altBase'],
    'ac ac cc':['altBase','refBase','altBase','refBase'],
    'ac ac ac':['E', 'E', 'E', 'E'],
    'ac ac aa':['refBase','altBase','refBase','altBase'],
    'ac aa ac':['altBase','refBase','refBase','refBase'],
    'ac aa aa':['refBase','altBase','refBase','refBase'],
    'aa cc ac':['refBase','refBase','altBase','altBase'],
    'aa ac ac':['refBase','refBase','altBase','refBase'],
    'aa ac aa':['refBase','refBase','refBase','altBase'],
    'aa aa aa':['refBase','refBase','refBase','refBase'],
    'aa a aa' :['refBase','refBase','refBase','refBase'],
    'ac c cc' :['altBase','refBase','altBase','altBase'],
    'ac a aa' :['refBase','altBase','refBase','refBase'],
    'cc a ac' :['altBase','altBase','refBase','refBase'],
    'cc c cc' :['altBase','altBase','altBase','altBase'],
    'aa c ac' :['refBase','refBase','altBase','altBase'],
    'ac c ac' :['refBase','altBase','altBase','altBase'],
    'ac a ac' :['altBase','refBase','refBase','refBase'],
    'aa a a'  :['refBase','refBase','refBase','refBase'],
    'cc a c'  :['altBase','altBase','refBase','refBase'],
    'cc c c'  :['altBase','altBase','altBase','altBase'],
    'ac a c'  :['altBase','refBase','refBase','refBase'],
    'ac c c'  :['altBase','refBase','altBase','altBase'],
    'ac c a'  :['refBase','altBase','altBase','altBase'],
    'ac a a'  :['refBase','altBase','refBase','refBase'],
    'aa c a'  :['refBase','refBase','altBase','altBase']
}

strange = {
'aa aa ac':'d',
'aa aa cc':'d',
'aa ac cc':'o',
'aa cc aa':'o',
'aa cc cc':'o',
'ac aa cc':'o',
'ac cc aa':'o',
'cc aa aa':'o',
'cc aa cc':'o',
'cc ac aa':'o',
'cc cc aa':'o',
'cc cc ac':'d'
}

def createHaploThreads(numP, numCh, GS, maxPos=len(snpMap)):
    res = [[] for i in range(numCh)]
    M = np.identity(5,dtype=np.int)
    baseMap = {b:i for i,b in enumerate(list('ACGT0'))}
    positions = []
    for k in range(maxPos):
        refA = snpMap[k][2] # majorA from filteredPos.txt
        altA = snpMap[k][3] # minorA from filteredPos.txt
        refM = {altA:'c', refA:'a', '0':'0'}
        rr = {'altBase':altA, 'refBase':refA}
        g = [ [GS[i,2*k], GS[i,2*k+1]] for i in range(numP)]

        # for example for trio
        # refA = 'C', altA = 'A'
        # refM = {'A': 'c', 'C': 'a', '0': '0'}
        # g=[['A','C'],['C','C'], ['C','C']]
        #

        counts = np.array([M[baseMap[x]] for i in range(2) for x in g[i]])

        # for example for g above
        # counts = array([[1, 0, 0, 0, 0],
        #                 [0, 1, 0, 0, 0],
        #                 [0, 1, 0, 0, 0],
        #                 [0, 1, 0, 0, 0]])
        #

        C = np.sum(counts, axis=0)

        # for example for counts above
        # C = array([1, 3, 0, 0, 0])
        # 
        positions.append(snpMap[k][:2] + [refA] + list( C[:4] ))
        for ch in range(numCh):
            gch = [g[0],g[1],g[2+ch]]
            if '0' in ''.join(g[0] + g[1] + g[2+ch]):
                res[ch].append(['c']*4)
                # 'c' for '0' genotype in any family member
                continue
            try:
                K = [sorted([refM[x[0]],refM[x[1]]]) for x in gch]
                # for example, K = [['a', 'c'], ['a', 'a'], ['a', 'a']]
                K = ' '.join([x[0] + x[1] for x in K])
                # for example, K = 'ac aa aa'
            except KeyError:
                print("k", k, "snpMap", ' '.join(snpMap[k]), "refM", refM, "gch", gch, file=sys.stderr)
                res[ch].append(['d']*4)
                # 'd' for denovo (or error in genotyping)
                continue
            if K in gn2ph:
                ht = gn2ph[K]
                if 'E' in ht:
                    b = ambiguous[gch[0][0]+gch[0][1]]
                    res[ch].append([b]*4)
                else:
                    ht = [rr[x] for x in ht]
                    res[ch].append(ht)
            else:
                # any violation of mendelian rule (either 'd' or 'o')
                res[ch].append([strange[K]]*4)
                          
    return positions, np.array(res)

positions, res = createHaploThreads(len(family), len(family) - 2, GS)

res1 = np.array(res).transpose()

HDB=outDir + '/'+famId

with open(HDB +'-pos.txt', 'w') as f:
    posHead='Chromosome Position RefAllele A C G T'.split(' ')
    f.write('\t'.join(posHead) + '\n')
    for p in positions:
        f.write('\t'.join(map(str,p))+'\n')
                
with open(HDB+'-hpth.txt', 'w') as f:
        personIds = members
        f.write( '\t'.join('Family childId haplothreadId Haplothread'.split(' ')) + '\n')
        htIds = 'MT MNT FT FNT'.split(' ')
        s = res1.shape
        for n in range(s[2]):
            for k in range(s[0]):
                f.write( '\t'.join([famId, personIds[2+n], htIds[k]] + [''.join(res1[k,:,n])]) + '\n')

