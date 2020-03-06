#!/usr/bin/env python
import numpy as np
from collections import defaultdict
import sys, os
import resource
import h5py

from GenomeAccess import openRef
if len(sys.argv) < 3:
    print "Usage: chip2HT.py <ped file> <bim file> [pedMap]" 
    sys.exit(1)

""" 

We process chunk of ped file and output chunk of ped file 
with sorted positions and eliminated bad positions
at the same time we write *-hpth.txt file for each family separately
and *-pos.h5 file for each family separately.
This is important, since we later aggregate families according to type:
TRIOS, QUADS, and MULTIPLEX
chip-sorted.bim is written when processing the first family

"""

pedFn = sys.argv[1]
bimFn = sys.argv[2]
chunkN = pedFn.split("/")[1].split(".")[0]

print >>sys.stderr,  pedFn, bimFn, chunkN

pMap = False
#"""
if len(sys.argv) >3:
    personFn = sys.argv[3]
    pMap=True
#"""

snpMap = []
persons = {}
nucFams = defaultdict(list)
Gmap={'1':'M', '2':'F', '0':'0'}
S = {'mom':0,'dad':1, 'prb':2,'sib':3}
cur_fam = None
FAM = {}

pedMapH=tuple('pedSampleId familyId personId fatherId motherId gender affectedStatus'.split(' '))
if pMap:
    personMap=defaultdict()
    personMap['0']='0'
    with open(personFn, 'r') as f:
        firstLine=f.readline().strip('\n\r').split('\t')
        assert (tuple(firstLine) == pedMapH)
        for l in f:
            cs = l.strip('\n\r').split('\t')
            personMap[cs[0]]=cs[1:]

CompM = {x:y for x,y in zip(list('ACGT0ID'),list('TGCA0ID'))}

chrM = {'chr'+str(v):i+1 for i,v in enumerate(range(1,23))}
snpIdx = defaultdict(int)
with open(bimFn, 'r') as f:
    n = 0
    for l in f:
        cs = l.strip('\n\r').split("\t")
        ch = cs[0]
        p = cs[3]
        snpIdx[(ch,p)] +=1

        if cs[3] == '-1':
            snpIdx[(ch,p)] +=1        

        snpMap.append(cs)
        n +=1

mapIdx = sorted([[i,[chrM[snpMap[i][0]],int(snpMap[i][3])]]
                 for i,v in enumerate(snpMap) if snpIdx[(snpMap[i][0],snpMap[i][3])] == 1], key=lambda x: x[1])

bimOutFn = 'chip-sorted.bim'
if not os.path.isfile(bimOutFn):
    outbim = open(bimOutFn, 'w')
    def writeOutBim(outbim):
        for m in mapIdx:
            k=m[0]    
            outbim.write('\t'.join(map(str,snpMap[k])) +'\n')
        outbim.close()

    writeOutBim(outbim)


chip_file='NF-HDB/'+chunkN + '.ped'
outchip = open(chip_file, 'w')

ambiguous = {
    'AC':'M', 'CA':'M', 
    'AG':'R', 'GA':'R', 
    'AT':'W', 'TA':'W', 
    'CG':'S', 'GC':'S', 
    'CT':'Y', 'TC':'Y',
    'GT':'K', 'TG':'K'
}

#    'mom dad child':[MT, MNT, FT, FNT]
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
#    'mom dad child':[MT, MNT, FT, FNT]
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

positions = defaultdict(list)
def createHaploThreads(numP, numCh, GS):
    #print >>sys.stderr, "maxPos", maxPos
    res = [[] for i in range(numCh)]
    M = np.identity(5,dtype=np.int)
    baseMap = {b:i for i,b in enumerate(list('ACGT0'))}
    #print >>sys.stderr, "id chr pos refA altA minorA majorA"
    for m in mapIdx:
        k=m[0]
        #print >>sys.stderr, k
        ch = snpMap[k][0]
        p = snpMap[k][3]

        # if we use real refA and do not do flip, we get tons of denovos
        #refA = snpMap[k][-1]
        
        # since we do not make flip, we use majorA as refA
        refA = snpMap[k][5]
        if snpMap[k][4] == '0':
            altA = refA 
        else:
            altA = snpMap[k][4]
            
        refM = {altA:'c', refA:'a', '0':'0'}
        rr = {'altBase':altA, 'refBase':refA}

        g = [ [GS[i,2*k], GS[i,2*k+1]] for i in range(numP)]

        counts = np.array([M[baseMap[x]] for i in range(2) for x in g[i]])
        C = np.sum(counts, axis=0)
        positions[k]=[snpMap[k][0], snpMap[k][3], refA] + list(C[:4])
        
        for ch in range(numCh):
            # g[0] - dad, g[1] - mom, g[2+ch] - child
            gch = [g[0],g[1],g[2+ch]]
            if '0' in ''.join(g[0] + g[1] + g[2+ch]):
                res[ch].append(['c']*4)
                continue
            try:
                K = [sorted([refM[x[0]],refM[x[1]]]) for x in gch]
                K = ' '.join([x[0] + x[1] for x in K])
            except KeyError:
                print >>sys.stderr, "k", k, "snpMap", ' '.join(snpMap[k]), "refM", refM, "gch", gch, "refA", refA, "altA", altA
                res[ch].append(['d']*4)
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
                res[ch].append([strange[K]]*4)
    return np.array(res)

class P():
    pass

def processFam(fams):
    for fid, fam in fams.items():
        mom = None
        dad = None
        if len(fam) <3:
            return False
        for v in fam:
            p = P();
            p.mom = v[3]
            if p.mom != '0':
                mom = p.mom
                p.role = 'prb' if v[5] == '2' else 'sib'                
            p.dad = v[2]
            if p.dad != '0':
                dad = p.dad
                p.role = 'prb' if v[5] == '2' else 'sib'
            p.gender = Gmap[v[4]]
            p.personId = v[1]
            p.famId = v[0]
            p.aff = v[5]
            nucFams[fid].append(p)
            persons[p.personId] = p
        if not mom or not dad:
            return False
        for p in nucFams[fid]:
            if p.personId == mom:
                p.role = 'mom'
            if p.personId == dad:
                p.role = 'dad'
        nucFams[fid] = sorted(nucFams[fid], key=lambda x: S[x.role])
    return True

def printPOS(cur_fam, positions, mapIdx):
    pos_file='NF-HDB/'+cur_fam+'-pos.h5'
    data_chr = np.array([positions[ip[0]][0] for ip in mapIdx])
    data_pos = np.array([positions[ip[0]][1] for ip in mapIdx], dtype=int)
    data_refA = np.array([positions[ip[0]][2] for ip in mapIdx])
    data_cnt = np.array([positions[ip[0]][3:] for ip in mapIdx], dtype=int)
    hf = h5py.File(pos_file, 'w')
    hf.create_dataset('chr', data=data_chr)
    hf.create_dataset('pos', data=data_pos)
    hf.create_dataset('cnt', data=data_cnt)
    hf.create_dataset('refA', data=data_refA)
    hf.close()
    
def printHDB(res1, members):
    htIds = 'FNT FT MNT MT'.split(' ')
    hdb_hpth = open('NF-HDB/'+cur_fam+'-hpth.txt', 'w')
    hdb_hpth.write( '\t'.join('Family childId haplothreadId Haplothread'.split(' ')) + '\n')
    s = res1.shape
    for n in range(s[2]):
        for k in range(s[0]):
            hdb_hpth.write( '\t'.join([cur_fam, members[2+n].personId, htIds[k]] + [''.join(res1[3-k,:,n])]) + '\n')
    hdb_hpth.close()
    
members = []
FAM = {}
global res
def createHT(FAM, cur_fam):
    global res
    print >>sys.stderr, "createing HT", cur_fam
    res = createHaploThreads(len(members), len(members) - 2, GS)
    res1 = np.array(res).transpose()
    printHDB(res1,members)
    printPOS(cur_fam, positions, mapIdx)
    return FAM.clear(), defaultdict(list)

n = 0
fams = defaultdict(list)
with open(pedFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        fId, pId, faId, moId, sex, aff = cs[:6]
        #print >>sys.stderr, fId, pId
        if cur_fam and fId != cur_fam:
            #print 'A'
            if processFam(fams):
                members = nucFams[cur_fam]
                for p in members:
                    FAM[p.personId]['role']=p.role
                GS = np.vstack(tuple([FAM[p.personId]['genotype'] for p in members]))
                createHT(FAM, cur_fam)
            else:
                print >>sys.stderr, "strange family", str(fams)
            FAM = {}
            fams = defaultdict(list)
            persons = {}
        #print 'B'
        cur_fam = fId
        if pMap:
            fId, pId, faId, moId, gender, aff = personMap[pId]
            cur_fam=fId
        fams[fId].append([fId, pId, faId, moId, sex, aff])
        G = cs[6:]
        Gfixed = [G[2*m[0]+i] for m in mapIdx for i in range(2)]
        outchip.write('\t'.join(cs[:6] + Gfixed) + '\n')
        FAM[pId] = {'personId':pId,'genotype':np.array(G)}
        n +=1
        #if n == 62:
        #    break

# process the last family
if processFam(fams):
    members = nucFams[cur_fam]
    for p in members:
        FAM[p.personId]['role']=p.role
    GS = np.vstack(tuple([FAM[p.personId]['genotype'] for p in members]))
    createHT(FAM, cur_fam)
else:
    print >>sys.stderr, "strange family", str(fams)

outchip.close()
print >>sys.stderr, "Done"

