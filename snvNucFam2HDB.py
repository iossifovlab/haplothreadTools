#!/bin/env python
from math import *
import numpy as np
import glob 
from time import time
import sys
import os.path
import pysam,gzip
import scipy.stats as sci

if len(sys.argv) < 3:
    print "Usage: snvNucFam2HDB.py <snvNucFamDir> <new HDB>  [chromosome]"
    sys.exit(1)
snvNucFamDir,HDB = sys.argv[1:3]

ch = None
if len(sys.argv) > 3:
   ch = sys.argv[3] 

print "snvNucFamDir:", snvNucFamDir, "HDB:", HDB, "ch:", ch

# genStatsF = sys.stdin

if os.path.isfile(snvNucFamDir + "/quadGenStat.bgz") and os.path.isfile(snvNucFamDir + "/baseCntStat.bgz"):
    if ch:
        GF = pysam.TabixFile(snvNucFamDir + "/quadGenStat.bgz")
        genStatsF = GF.fetch(region="%s:1-500000000" % ch)
        BF = pysam.TabixFile(snvNucFamDir + "/baseCntStat.bgz")
        baseCntStatsF = BF.fetch(region="%s:1-500000000" % ch)
    else:
        genStatsF = gzip.open(snvNucFamDir + "/quadGenStat.bgz") 
        baseCntStatsF = gzip.open(snvNucFamDir + "/baseCntStat.bgz") 
elif os.path.isfile(snvNucFamDir + "/quadGenStat.txt") and os.path.isfile(snvNucFamDir + "/baseCntStat.txt"): # this is for testing purposes
    genStatsF = open(snvNucFamDir + "/quadGenStat.txt")
    baseCntStatsF = open(snvNucFamDir + "/baseCntStat.txt")
else:
    print "No quadGenStat file found"
    sys.exit(1)

F = open(snvNucFamDir + '/mainNucFamDir/params.txt', 'r')
params = {}
for line in F:
    if line.startswith("#"):
        continue
    k,v = line.strip("\n\r").split("=")
    params[k] = v
F.close()


personIds = params["personIds"].split(',')
familyId = params["quad.quad_id"]
if os.path.isfile(snvNucFamDir + "/../../snvNucFamAnnot/" + familyId + "/snvDenovoSW.txt"):
    SDF = open(snvNucFamDir + "/../../snvNucFamAnnot/" + familyId + "/snvDenovoSW.txt", 'r')
    strongDN = {}
    for line in SDF:
        loc, bestState = line.strip('\n\r').split('\t')        
        cs = loc.strip('\n\r').split(':')        
        if not cs[0] == ch:
            continue
        elif not cs[1] in strongDN:
            strongDN[cs[1]] = map(int,bestState.split('/')[1].split(' '))
    SDF.close()
else:
    print "not file " + snvNucFamDir + "/../../snvNucFamAnnot/" + familyId + "/snvDenovoSW.txt found"
    exit(1)

#for k in strongDN:
#    print k, strongDN[k]

print "familyId:", familyId, "personIds:", personIds

numP = len(personIds)
numCh = numP-2
flatten = lambda l: [item for sublist in l for item in sublist]
#chIds = flatten([[k]*4 for k in personIds[2:]])


bsStateScoreMin = 50
chi2PvalMin = -4
maxAltFreqPrcntToCache = 5.0


st = time()

bsList3 = [
'000/222',
'010/212',
'011/211',
'021/201',
'100/122',
'101/121',
'110/112',
'111/111',
'112/110',
'121/101',
'122/100',
'201/021',
'211/011',
'212/010',
'222/000'
]

c2g = {'02':'cc', '20':'aa', '11':'ac'}

bs2gn3 = {}

def createBs2gn(BS, BS2GN):
    for bs in BS:
        L = len(bs)/2
        BS2GN[bs] = ' '.join([c2g[bs[i]+bs[i+L+1]] for i in range(L)])
    
createBs2gn(bsList3,bs2gn3)

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
'''
matchMap ={(a+b):int(a==b) for a in list('ACGT') for b in list('ACGT')}
for b in list('ACGTED'):
    for a in list('ED'):
        matchMap[(b+a)] = 1
        matchMap[(a+b)] = 1
'''
# quadGenStat.bgz record structure
# 22      16050159        C,24,T,1121/1101,118,5,

if ch == None:
    specialPos = open(snvNucFamDir + "/log/specialPos.txt", 'w')
else:
    specialPos = open(snvNucFamDir + "/log/specialPos-" + ch +".txt", 'w')

M = {0:'A',1:'C',2:'G', 3:'T'}
RM = {v:k for k,v in M.items()}

'''
def processSpecial(chr, pos, refA, depth, altA, bSt, numP, numCh, cs, csCnt,e=0.01):
    GEN = {0:refA + refA, 1:refA+altA, 2:altA+altA}

    cnts = [map(float, k.split(',')) for k in csCnt[2:]]
    altCnt = cnts[0][RM[altA]] + cnts[1][RM[altA]] 
    refCnt = cnts[0][RM[refA]] + cnts[1][RM[refA]] 
    if (not refA == altA) and  \
       ((cnts[0][RM[refA]] + cnts[0][RM[altA]]) < 10 or \
        (cnts[1][RM[refA]] + cnts[1][RM[altA]]) < 10):
        #print cnts[0][RM[refA]], cnts[0][RM[altA]], cnts[1][RM[refA]], cnts[1][RM[altA]], 'low coverage'
        return [['c']*4]*numCh
    elif (refA == altA) and ((cnts[0][RM[refA]] < 10) or (cnts[0][RM[refA]] < 10)):
        return [['c']*4]*numCh
    ap = [1 - altCnt/(refCnt + altCnt), altCnt/(refCnt + altCnt)]
    gp = [ap[0]**2, 2*ap[0]*ap[1], ap[1]**2]
    E = [1-e, gp[1], 1-e]
    all = [refA, altA, altA]
    gens = []
    odds = []
    for l in [0,1]: 
        pH = [sci.binom.pmf(cnts[l][RM[all[k]]], cnts[l][RM[refA]] + cnts[l][RM[altA]], E[k])*gp[k] for k in range(3)]
        res = sorted([[v,k] for k,v in enumerate(pH)],reverse=True)
        #print refA, altA, bSt, cnts[:2], GEN[res[0][1]], res[1][0]/res[0][0], res
        gens.append(GEN[res[0][1]])
        odds.append(max(0.01, res[1][0]/res[0][0]))
    for chId in range(numCh):
        id = 2 + chId
        if (not refA == altA) and  \
           (cnts[id][RM[refA]] + cnts[id][RM[altA]]) < 10:
            return ['c']*4
        elif (refA == altA) and (cnts[id][RM[refA]] < 10):
            return ['c']*4
        bStTrio = bSt[:2] + bSt[id] + '/' + bSt[(numP + 1):(numP+3)] + bSt[numP + 1 + id]
        pH = [sci.binom.pmf(cnts[id][RM[all[k]]], cnts[id][RM[refA]] + cnts[id][RM[altA]], E[k])*gp[k] for k in range(3)]
        res = sorted([[v,k] for k,v in enumerate(pH)],reverse=True)
        gens.append(GEN[res[0][1]])
        odds.append(max(0.01,res[1][0]/res[0][0]))
        denovo = False
        omission = False
        if gens[0] == gens[1] and \
           (gens[0][0] == gens[0][1]) and \
           (not gens[id][0] == gens[id][1]):
            # AA,AA, AC
            denovo = True
        if (not gens[0] == gens[1]) and \
             (gens[0][0] == gens[0][1]) and \
             (gens[1][0] == gens[1][1]) and \
             (gens[id][0] == gens[id][1]):
            # AA, CC, CC
            # AA, CC, AA
            omission = True
        if denovo:
            gens.append('denovo')
        if omission:
            gens.append('omission')

    specialPos.write('\t'.join(cs + csCnt[2:] + [','.join(gens)] + [','.join(map(str,odds))]) + '\n')

def genotype(cnt):
    amap = {0:'A',1:'C',2:'G', 3:'T'}
    D = sum(cnt)
    if D < 20:
        return 'low'
    C = sorted([[v,i] for i,v in enumerate(cnt)], reverse=True)
    M = C[0][0]
    m = C[1][0]
    o = C[2][0] + C[3][0]
    Mf = M/D
    mf = m/D
    of = o/D
    if of > 0.01 and o > 1:
        return 'conf'
    elif mf > 0.05 or m < 3:
        return amap[C[0][1]]+amap[C[0][1]]
    elif mf > 0.25 and m > 5:
        a = sorted(amap[C[0][1]]+amap[C[1][1]])
        return a[0] + a[1]
    else:
        return 'conf'
'''

def genotypeTrio(cnt,refA,):
    amap = {0:'A',1:'C',2:'G', 3:'T'}
    RM = {v:k for k,v in amap.items()}
    cnt = [map(float, k) for k in cnt]
    D = [sum(k) for k in cnt]
    if min(D) < 15:
        return 'low'
    C = sorted([[v,i] for i,v in enumerate(cnt[2])], reverse=True)
    M = C[0][0]
    m = C[1][0]
    mf = m / D[2]
    mC = sorted([[v,i] for i,v in enumerate(cnt[0])], reverse=True)
    fC = sorted([[v,i] for i,v in enumerate(cnt[1])], reverse=True)
    mA = mC[0][0] == D[0]
    fA = fC[0][0] == D[1]
    Mhom = cnt[0][RM[refA]] == D[0]
    Fhom = cnt[1][RM[refA]] == D[1]
    #if Mhom and Fhom and ((not (amap[C[1][1]] == refA)) and (mf > 0.15 and m >6) or (not (amap[C[0][1]] == refA))) :
    #    return 'denovo'
    if (mA and cnt[2][mC[0][1]] == 0) or (fA and cnt[2][fC[0][1]] == 0):
        return 'omission'
    else:
        return 'conf'

def processSpecial(chr, pos, refA, depth, altA, bSt, numP, numCh, cs, csCnt,e=0.01):
    cnts = [map(float, k.split(',')) for k in csCnt[2:]]
    mD = sum(cnts[0])
    fD = sum(cnts[1])
    gens = []
    for n in range(numCh):
        id = 2 + n
        gens.append( genotypeTrio([cnts[0], cnts[1], cnts[id]], refA) )
        specialPos.write('\t'.join(cs + csCnt[2:] + [gens[n]])  + '\n')
    return gens

def printHaploThreads(numP, numCh):
    allMap = {v:int(k) for v,k in zip(list('ACGT'),list('0123'))}
    res = [[] for i in range(numCh)]
    positions = []
    n = 0
    hd = 'chrom pos refA ,depth, altA, bestState, varScr, bsScr, chi2Scr'.split(' ')
    for line, cntLine in zip(genStatsF, baseCntStatsF):
        cs = line.strip().split('\t')
        csCnt = cntLine.strip().split('\t')
        chr,pos = cs[:2]
        assert(cs[:2]==csCnt[:2])
        refA, depth, altA, bSt, varScr, bsScr, chi2Scr = cs[2].split(',')
        if int(bsScr) < 50 or float(chi2Scr) < chi2PvalMin:
            gens = processSpecial(chr, pos, refA, depth, altA, bSt, numP, numCh, cs, csCnt)
            positions.append([chr, pos, refA] + [0,0,0,0])
            for chId in range(numCh):
                id = 2 + chId
                #print pos, type(pos), chId, type(chId)
                #print (pos in strongDN)
                #print (strongDN[pos][id] == 1)
                if pos in strongDN and strongDN[pos][id] == 1:
                    specialPos.write('\t'.join(cs + csCnt[2:] + ['AAA from snvDenovoSW'])  + '\n')
                    res[chId].append(['d']*4)
                    continue
                if 'low' in gens[chId]:
                    ht = ['c']*4
                    res[chId].append(ht)
                elif 'conf' in gens[chId]:
                    ht = ['u']*4
                    res[chId].append(ht)
                elif 'omission' in gens[chId]:
                    ht = ['o']*4
                    res[chId].append(ht)
                #elif 'denovo' in gens[chId]:
                #    #ht = ['d']*4
                #    ht = ['u']*4
                #    res[chId].append(ht)
                else:
                    print >>sys.stderr, 'Strange genotype ', gens[chId]
                n=n+1
                #if n >1000:
                #    break
            continue
        pRec = [chr, pos, refA]
        cnt = [0,0,0,0]
        cnt[allMap[altA]] = cnt[allMap[altA]] + int(bSt[numP + 1]) + int(bSt[numP + 2])
        cnt[allMap[refA]] = cnt[allMap[refA]] + int(bSt[0]) + int(bSt[1])
        positions.append( pRec + cnt )
        rec = {}
        rec['E'] = 'E'
        rec['refBase'] = refA
        rec['altBase'] = altA
        for chId in range(numCh):
            id = 2 + chId
            if pos in strongDN and strongDN[pos][id] == 1:
                specialPos.write('\t'.join(cs + csCnt[2:] + ['BBB from snvDenovoSW'])  + '\n')
                res[chId].append(['d']*4)
                continue
            bStTrio = bSt[:2] + bSt[id] + '/' + bSt[(numP + 1):(numP+3)] + bSt[numP + 1 + id]
            ht = [rec[name] for name in gn2ph[bs2gn3[bStTrio]]]
            if 'E' in ht:
                ht = [ambiguous[refA+altA]]*4
            res[chId].append(ht)
        n = n+1
    return res, positions

res, positions = printHaploThreads(numP, numCh)
specialPos.close()

if ch == None:
    FHT = open(HDB + '-hpth.txt', 'w')
else:
    FHT = open(HDB + ch +'-hpth.txt', 'w') 
htIds = 'MT MNT FT FNT'.split(' ')

FHT.write( '\t'.join('Family childId haplothreadId Haplothread'.split(' ')) + '\n')
res1 = np.array(res).transpose()
s = res1.shape
for n in range(s[2]):
    for k in range(s[0]):
        FHT.write( '\t'.join([familyId, personIds[2+n], htIds[k]] + [''.join(res1[k,:,n])]) + '\n')
FHT.close()
if ch == None:
    FPOS = open(HDB + '-pos.txt', 'w')
else:
    FPOS = open(HDB + ch +'-pos.txt', 'w')
FPOS.write( '\t'.join('Chromosome Position RefAllele A C G T'.split(' ')) + '\n')
for line in positions:
    FPOS.write( '\t'.join(map(str,line)) + '\n')
FPOS.close()



