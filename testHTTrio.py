#!/bin/env python
from math import *
import numpy as np
import sys
from metaData import *
from collections import Counter

if len(sys.argv) < 4:
    print "Usage: testHTTrio.py <htTrioFile> <pedTrioFile> <htTrioPosFn>"
    sys.exit(1) 

htTrioFn = sys.argv[1] 
pedTrioFn = sys.argv[2] 
htTrioPosFn = sys.argv[3] 

htTrio = {}
with open(htTrioFn, 'r') as f:
    f.readline()
    for l in f:
        cs = l.strip('\n\r').split('\t')
        htTrio[cs[2]] = np.array(list(cs[3]))

refAA = []
with open(htTrioPosFn, 'r') as f:
    f.readline()
    for l in f:
        cs = l.strip('\n\r').split('\t')
        refAA.append(cs[2])

pedTrio = {}
with open(pedTrioFn, 'r') as f:
    fam = ['mom', 'dad', 'kid']
    for k in range(3):
        cs = f.readline().strip('\n\r').split(' ')
        tmp = np.array(cs[6:]).reshape((len(cs[6:])/2,2))
        pedTrio[fam[k]] = [''.join(sorted(x)) for x in tmp]
pedTrio['GN'] = [x for x in zip(pedTrio['mom'], pedTrio['dad'], pedTrio['kid'])]

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
    'aa aa ac':['d','d','d','d'],
    'cc cc ac':['d','d','d','d']
}

ambiguous = {
    'AC':'M', 'CA':'M', 
    'AG':'R', 'GA':'R', 
    'AT':'W', 'TA':'W', 
    'CG':'S', 'GC':'S', 
    'CT':'Y', 'TC':'Y',
    'GT':'K', 'TG':'K'
}

def makeHT(refA, GN, th, n):
    """ refA is reference allele, one of list('ACGT')
    GN is a list of genotypes of mom, dad, and child, where genotype is 
    concatenation of two alleles
    th is one of ['MT, 'MNT', 'FT', 'FNT']
    n is a number
    """
    #print >>sys.stderr, "refA: ", refA
    #print >>sys.stderr, "GN: ", GN
    #print >>sys.stderr, "th: ", th
    #print >>sys.stderr, "n: ", n    
    gm = {'A': 'c', 'C': 'c', 'G': 'c', 'T': 'c'}
    htm = {'MT':0, 'MNT':1, 'FT':2, 'FNT':3}
    gm[refA] = 'a'
    if '0A' in GN or '0C' in GN or '0G' in GN or '0T' in GN  or '00' in GN:
        return 'c'
    if len(set(GN)) == 1:
        if len(set(list(GN[0]))) == 1:
            return GN[0][0]
    alleles = set(list(''.join(GN)))
    #print 'alleles: ', alleles
    altA = [i for i in alleles if i != refA]
    if len(altA) >1:
        return 'u'
    #altA = altA[0]    
    altA = altA[0]
    #print 'altA: ', altA
    mHom = len(set(list(GN[0]))) == 1
    dHom = len(set(list(GN[1]))) == 1
    #print 'mHom, dHom: ', mHom, dHom
    if mHom and (not GN[0][0] in GN[2]) or dHom and (not GN[1][0] in GN[2]):
        return 'o'
    if (not GN[2][0] in GN[0]) and (not GN[2][1] in GN[0]) or \
       (not GN[2][0] in GN[1]) and (not GN[2][1] in GN[1]):
        return 'd'

    gg = ['']*3
    for n in range(3):
        g = sorted([gm[i] for i in GN[n]])
        gg[n] = g[0] + g[1]
    
    r = gn2ph[' '.join(gg)][htm[th]]
    #print 'GN, gg, th, r; ', GN, gg, th, gn2ph[' '.join(gg)]
    if r == 'refBase':
        return refA
    elif r == 'altBase':
        return altA
    elif r == 'E':
        return ambiguous[GN[0]]
    elif r=='d':
        return 'd'
                
for h in 'MT MNT FT FNT'.split(' '):
    pedTrio[h] = np.array([makeHT(refAA[k], pedTrio['GN'][k], h, 1) for k in range(len(pedTrio['GN']))])

compare = {}
for h in 'MT MNT FT FNT'.split(' '):
    compare[h] = np.array([int(x[0] == x[1]) for x in zip(htTrio[h], pedTrio[h])] )
    print 100.*sum(compare[h])/len(compare[h])

idx = np.where(compare['MT'] == 0)[0]

print Counter(htTrio['MT'][idx])
print Counter(pedTrio['MT'][idx])

a = [[htTrio['MT'][k], pedTrio['MT'][k], pedTrio['GN'][k]] for k in idx]
for x in a:
    print ' '.join(x[:2] + list(x[2]))
