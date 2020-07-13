#!/usr/bin/env python
import os, sys
from collections import defaultdict, Counter
import pickle
from glob import glob
from time import time
import numpy as np
import scipy.stats as sc
from hw import *
from GenomeAccess import openRef
import json

with open("config.json") as f:
    config=json.load(f)

par1par2=config["parameters"]["par1"]+"," + config["parameters"]["par2"]


"""
chr  
pos  
refA  
majorA  
minorA  
MM_N  
Mm_N  
mm_N 
00_N  
majorA_N = (2*MM_N + Mm_N) 
minorA_N = (2*mm_N + Mm_N) 
0_N      = (2*OO_N) 
percentCalled = 100*(majorA_N + minorA_N)/(majorA_N + minorA_N + 0_N)  
MAF (minor allele frequency) = minorA_N / (majorA_N + minorA_N) 
HW = hw_test(MM_N,Mm_N,mm_N) 
"""

if len(sys.argv) < 5:
    print("Usage: chipStatSummary.py  <statSummary.stat file> <chip.map file>> <number of families> [<agre pos file>]")
    sys.exit(1)

summaryFn = sys.argv[1]
mapFn = sys.argv[2]
sscPosFn = sys.argv[3]
famN=int(sys.argv[4])

print(summaryFn, file=sys.stderr)
print(mapFn, file=sys.stderr)
print(sscPosFn, file=sys.stderr)
print(famN, file=sys.stderr)

agre = False
if len(sys.argv) >5:
    agrePosFn = sys.argv[5]
    print(agrePosFn, file=sys.stderr)
    agre = True
    
with open(summaryFn, 'rb') as F:
    stats = np.load(F,allow_pickle=True)

N = famN*2
parents = N

snpMap = []

compM = {x:y for x,y in zip(list('ACGT0ID'),list('TGCA0ID'))}
chroms = ['chr' + str(i) for i in range(1,24)]
snpDict = defaultdict(int)

st = time()
with open(mapFn, 'r') as f:
    n = 0
    for l in f:
        cs = l.strip('\n\r').split('\t')
        ch = 'chr' +cs[0]
        p = cs[3]
        snpMap.append(cs)
        snpDict[(ch,p)] += 1
        n +=1
print("Read chip.map", time() - st, file=sys.stderr)

sscPos = {}
bases=list('ACGT0')
with open(sscPosFn, 'r') as f:
    for l in f:
        if l.startswith('C'):
            continue
        cs = l.strip('\n\r').split('\t')
        if not (cs[0],cs[1]) in snpDict:
            continue
        cnt = sorted([[bases[i-3], int(cs[i])] for i in range(3,7) 
                      if int(cs[i]) > 0], key=lambda x: x[1], 
                     reverse=True )

        sscPos[(cs[0],cs[1])] = ''.join([x[0] for x in cnt])

print("Read ssc positions", time() - st, file=sys.stderr)

if agre:
    agrePos = {}
    with open(agrePosFn, 'r') as f:
        for l in f:
            if l.startswith('C'):
                continue
            cs = l.strip('\n\r').split('\t')
            if not (cs[0],cs[1]) in snpDict:
                continue
            cnt = sorted([[bases[i-3], int(cs[i])]
                          for i in range(3,7) 
                          if int(cs[i]) > 0], key=lambda x: x[1], 
                     reverse=True )

            agrePos[(cs[0],cs[1])] = ''.join([x[0] for x in cnt])

    print("Read agre positions", time() - st, file=sys.stderr)

MAP = [[x[0], x[1], x[3], x[4]] for x in snpMap]
mapLen = len(MAP)


def verify(sP, k):
    basesP = [item for sublist in sP for item in sublist]
    c = {x:sum([sP[y]*2 if (x,x) == y else sP[y]
                for y in sP if x in y]) for x in basesP}
    return sum(c.values())

def alleles(sPF, sPM, sC, k):
    sP =  sPF + sPM
    basesP = [item for sublist in sP for item in sublist]
    c = {x:sum([sP[y]*2 if (x,x) == y else sP[y]
                for y in sP if x in y]) for x in basesP}
    d = sorted(list(c.items()), key=lambda x: x[1], reverse=True)
    d = [x for x in d if x[0] != '0']
    if len(d) ==0:
        return ['0', '0', 0, 0, 0, N, 0, 0, N, 0, 0, 0, N, 0, 0, 0, N, 0,0,N, 0,0,0]

    majorA = d[0][0]
    
    if len(d)>1:
        minorA = d[1][0]
    else:
        minorA = majorA
        # 
        # look if  children have another alleles
        basesC = [item for sublist in sC for item in sublist]        
        c = {x:sum([sC[y]*2 if (x,x) == y else sC[y] for y in sC if x in y])
             for x in basesC}
        d = sorted(list(c.items()), key=lambda x: x[1], reverse=True)
        #print >>sys.stderr, '\t'.join(map(str,[MAP[k][0],MAP[k][1],
        #                                 'parents', sP, 'children', sC, d]))
        d = [x for x in d if x[0] != '0']
        
        if len(d) > 1:
            minorA = d[1][0]
        if len(d) > 2:
            print("Strange", d, file=sys.stderr)
    
    keys = [[(majorA, majorA)],
            [(minorA, minorA)],
            [(majorA, minorA), (minorA, majorA)],
            [('0','0')]]
    ### *_NF  for dad and *_NM for mom
    MM_NF, mm_NF, Mm_NF, OO_NF = [sum([sPF[i] for i in k if i in sPF]) for k in keys]
    #print >>sys.stderr, majorA, minorA, MM_NF, mm_NF, Mm_NF, OO_NF
    MM_NM, mm_NM, Mm_NM, OO_NM = [sum([sPM[i] for i in k if i in sPM]) for k in keys]
    #print >>sys.stderr, majorA, minorA, MM_NM, mm_NM, Mm_NM, OO_NM
    MM_NC, mm_NC, Mm_NC, OO_NC = [sum([sC[i] for i in k if i in sC]) for k in keys]    
    OO_N = OO_NF + OO_NM
    MM_N = MM_NF + MM_NM
    Mm_N = Mm_NF + Mm_NM
    mm_N = mm_NF + mm_NM

    MM_NT = MM_N + MM_NC
    mm_NT = mm_N + mm_NC
    Mm_NT = Mm_N + Mm_NC
    OO_NT = OO_N + OO_NC

    if minorA == majorA:
        mm_N = 0
        Mm_N = 0
    majorA_N = 2*MM_N + Mm_N
    minorA_N = 2*mm_N + Mm_N
    O_N = 2*OO_N
    majorA_NT = 2*MM_NT + Mm_NT
    minorA_NT = 2*mm_NT + Mm_NT
    O_NT = 2*OO_NT

    percentCalled = '%.2f' % (100.0*(majorA_NT + minorA_NT)/(majorA_NT + minorA_NT + O_NT)) if majorA_NT + minorA_NT + O_NT > 0 else '0.00'


    if majorA == minorA  or sum(HW1(MM_N, Mm_N,mm_N)==0):
        p = 1.0
    else:
        """
        if (MAP[k][0] == '23' and not isPseudoAutosomalX( MAP[k][2], par1par2 )):
            #print (pos, "A", file=sys.stderr)
            hw, p = sc.chisquare(np.array([MM_NF, mm_NF, MM_NM, Mm_NM,mm_NM]),
                                   f_exp=HW_X(MM_NF, mm_NF, MM_NM, Mm_NM,mm_NM), ddof=3)
            print(MAP[k][2], MM_NF, mm_NF, MM_NM, Mm_NM,mm_NM, hw, p, file=sys.stderr)
        else:
            #print (, "B", file=sys.stderr)
            hw, p = sc.chisquare(np.array([MM_N, Mm_N,mm_N]), f_exp=HW1(MM_N, Mm_N,mm_N), ddof=1)
        """
        if ((MAP[k][0] == '23' or MAP[k][0] == 'chrX') and not isPseudoAutosomalX( MAP[k][2], par1par2 )):
            #print (pos, "A", file=sys.stderr)
            p, eCnt = pval_count_X([MM_NF, mm_NF, MM_NM, Mm_NM,mm_NM])
            #print(MAP[k][2], MM_NF, mm_NF, MM_NM, Mm_NM,mm_NM, hw, p, file=sys.stderr)
        else:
            #print (, "B", file=sys.stderr)
            p, eCnt = pval_count_autosome([MM_N, Mm_N,mm_N])

    p = '%.5f' % p

    indel_N = N*4 - (majorA_N + minorA_N + O_N)

    if majorA != minorA:
        MAF = '%.4f' % (1.0*minorA_N / (majorA_N + minorA_N))
        return [majorA, minorA, MM_NT, Mm_NT, mm_NT, OO_NT,  majorA_NT, minorA_NT, O_NT, MM_NF, Mm_NF, mm_NF, OO_NF, MM_NM, Mm_NM, mm_NM, OO_NM, majorA_N, minorA_N, O_N, percentCalled,MAF,p]
    else:
        MAF = 0.0
    return [majorA, '', MM_NT, Mm_NT, mm_NT, OO_NT,  majorA_NT, minorA_NT, O_NT, MM_NF, Mm_NF, mm_NF, OO_NF, MM_NF, Mm_NF, mm_NF, OO_NF,  majorA_N, minorA_N, O_N, percentCalled,MAF,p]

statsPF = stats[0]
statsPM = stats[1]
statsC  = stats[2]

for k in range(statsPF.shape[0]):
    CNT = verify(statsPF[k]+statsPM[k], k)
    if  CNT !=  N*2:
        print(k, statsPF[k]+statsPM[k], CNT, "!=", N*2, parents, file=sys.stderr)

def rejection(sscAls, inSSC, ALS, ch, pos, agreAls="", inAGRE=False, agre=False):

    sscA = set(list(sscAls))
    sscAcomp = set([compM[x] for x in sscA])
    chipA = set(list(ALS['majorA']+ALS['minorA']))
    chipAcomp = set([compM[x] for x in chipA ])
    indel = ('I' in chipA or 'D' in chipA or
             len(ALS['majorA']) > 1 or len(ALS['minorA']) >1)
    sscBiallelic = len(sscAls) == 2
    if agre:
        agreA = set(list(agreAls))
        agreAcomp = set([compM[x] for x in agreA])
        agreBiallelic = len(agreAls) == 2
        
    flip = 'yes' if (inSSC and sscBiallelic 
                     and chipAcomp <= sscA 
                     and sscA != sscAcomp and not indel) else 'no'
    
    if flip == 'yes' and inSSC and sscBiallelic and not indel:
        print("flip", ch, pos, "in chip", chipA, "in ssc", sscA, file=sys.stderr)

    if agre:
        compatible = True if (inSSC 
                            and sscBiallelic
                            and inAGRE
                            and agreBiallelic
                            and ((sscA <= chipA or sscAcomp <= chipA) 
                            or (chipA <= sscA or chipAcomp <= sscA))
                            and ((agreA <= chipA or agreAcomp <= chipA) 
                            or (chipA <= agreA or chipAcomp <= agreA))
                  and not indel) else False
    else:
        compatible = True if (inSSC 
                  and sscBiallelic
                  and ((sscA <= chipA or sscAcomp <= chipA) 
                       or (chipA <= sscA or chipAcomp <= sscA))
                  and not indel) else False

    if agre:
        if not compatible and not indel and sscBiallelic and  inSSC and inAGRE and agreBiallelic:
            print("not compatible", ch, pos, "in chip", chipA, "in ssc", sscA, "in agre", agreA, file=sys.stderr)
    else:
        if not compatible and not indel and sscBiallelic and  inSSC:
            print("not compatible", ch, pos, "in chip", chipA, "in ssc", sscA, file=sys.stderr)
        
    reason = ""
    if indel:
        reason += "indel;"
    if snpDict[(ch, pos)] > 1:
        reason += "multiple pos;"

    if float(ALS['HW_p']) < 0.00001:
        reason += "HW < 0.00001;"
    
    if float(ALS['percentCalled']) < 95:
        reason += "percentCalled < 95;"
    if not inSSC:
        reason += "not in SSC;"
    if agre and not inAGRE:
        reason += "not in AGRE;"        
    if not sscBiallelic:
        reason += "ssc not biallelic;"
    if agre and not agreBiallelic:
        reason += "agre not biallelic;"        
    if float(ALS['MAF']) == 0:
        reason += "MAF is 0;"
    if not ch in chroms:
        reason += "not autosome;"
    if not compatible:
        reason += "not compatible;"

    return [flip, 'yes' if reason != '' else 'no', reason]

out = []
st = time()
if agre:
    hdAgre='chrom varId pos refA sscAlleles agreAlleles majorA minorA MM_N Mm_N mm_N 00_N majorA_N minorA_N 0_N MM_dad Mm_dad mm_dad OO_dad MM_mom Mm_mom mm_mom OO_mom majorA_parents minorA_parents O_parents percentCalled MAF HW_p inSSC inAGRE flip reject reason'.split(' ')
    print('\t'.join(hdAgre))
else:
    hd='chrom varId pos refA sscAlleles majorA minorA MM_N Mm_N mm_N 00_N majorA_N minorA_N 0_N MM_dad Mm_dad mm_dad OO_dad MM_mom Mm_mom mm_mom OO_mom majorA_parents minorA_parents O_parents percentCalled MAF HW_p inSSC flip reject reason'.split(' ')
    print('\t'.join(hd))


for k in range(statsPF.shape[0]):
    out = MAP[k]
    ch = 'chr'+MAP[k][0]
    p = MAP[k][2]
    present = sscPos[(ch, p)] if (ch,p) in sscPos else ""
    out += [present]
    if agre:
        presentA = agrePos[(ch, p)] if (ch,p) in agrePos else ""
        out += [presentA]
    als = alleles(statsPF[k],statsPM[k],statsC[k], k)
    if agre:
        ALS = {k:v for k,v in zip(hdAgre[6:-5], als)}
    else:
        ALS = {k:v for k,v in zip(hd[5:-4], als)}        
    out += als
    inSSC = ('chr'+MAP[k][0], MAP[k][2]) in sscPos
    out += [inSSC]
    #if agre:
    inAGRE = agre and ('chr'+MAP[k][0], MAP[k][2]) in agrePos
    out += [inAGRE]
    if agre:
        out += rejection(present, inSSC, ALS, ch, p,
                         presentA, inAGRE, agre)
    else:
        out += rejection(present, inSSC, ALS, ch, p)
    print('\t'.join(map(str,out)))
