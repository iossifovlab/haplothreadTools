#!/usr/bin/env python
import sys,os
import numpy as np
from collections import defaultdict, Counter 
import pickle


if len(sys.argv) <4:
    print "Usage: ht2chip.py <hdb> <original bim> <fms.pckl>"
    exit(1)

hdb = sys.argv[1]
obimFn = sys.argv[2]
fmsFn = sys.argv[3]

htFn = hdb + '-hpth.txt'
posFn = hdb + '-pos.txt'
pedFn = hdb + '.ped'
bimFn = hdb + '.bim'
famFn = hdb + '.fam'

ambiguous = {
    'AC':'M', 'CA':'M', 
    'AG':'R', 'GA':'R', 
    'AT':'W', 'TA':'W', 
    'CG':'S', 'GC':'S', 
    'CT':'Y', 'TC':'Y',
    'GT':'K', 'TG':'K'
}

amb = {v:list(k) for k,v in ambiguous.items()}
amb['d'] = ['0','0']
amb['c'] = ['0','0']
amb['o'] = ['0','0']
amb['u'] = ['0','0']

amb_k = amb.keys()

persons = defaultdict()

F = open(fmsFn)
families,badFamilies = pickle.load(F)
F.close()

genM = {'M':'1', 'F':'2'}

for fd in families.values():
    for pd in fd.memberInOrder:
        fId = fd.familyId
        pId = pd.personId
        dad = fId + '.fa' if pd.role in ['sib', 'prb'] else '0'
        mom = fId + '.mo' if pd.role in ['sib', 'prb'] else '0'
        persons[pId]=[fId, pId, dad, mom, genM[pd.gender], '2' if pd.role == 'prb' else '1'] 

OBIM = defaultdict()
with open(obimFn, 'r') as f:
    for l in f:
        cs = l.strip("\n\r").split("\t")
        OBIM[(cs[0],cs[3])] = cs

start = True
ALS = None
pedOut = open(pedFn, 'w')

def processFamily(family):
    global ALS
    global start
    print >>sys.stderr, "processing family", family[0]['fId']
    assert len(set([x['fId'] for x in family]))==1
    fId = family[0]['fId']
    assert set([len(x) for x in family]) == set([4])
    L = len(family[0]['ht'])
    #print >>sys.stderr, "ht length", L
    assert set([len(x['ht']) for x in family]) == set([L])
    mom = persons[fId + '.mo']
    dad = persons[fId + '.fa']

    children = sorted(set([x['pId'] for x in family 
                           if 'p' in x['pId'] 
                           or 's' in x['pId']]))

    #print >>sys.stderr, "children", children
    childG = np.empty((len(children), 2*L), dtype='S1')
    if start:
        ALS = np.array([Counter() for k in range(L)])
    for i,ch in enumerate(children):
        H = {x['htId']:x['ht'] for x in family if x['pId'] == ch}
        #print >>sys.stderr, "A"
        childG[i] = np.array([[H['FT'][k],H['MT'][k]] 
                              if not H['FT'][k] in amb_k
                                      else amb[H['FT'][k]]
                              for k in range(L)]).reshape((2*L))
        if i == 0:
            tmp = [[H['MT'][k],H['MNT'][k]]
                             if not H['MT'][k] in amb_k
                                      else amb[H['MT'][k]]
                             for k in range(L)]
            ALS +=  np.array([Counter(x) for x in tmp])
            momG = np.array(tmp).reshape((2*L))
            tmp = [[H['FT'][k],H['FNT'][k]]
                             if not H['FT'][k] in amb_k
                                      else amb[H['FT'][k]] 
                             for k in range(L)]
            ALS +=  np.array([Counter(x) for x in tmp])
            dadG = np.array(tmp).reshape((2*L))
            #print >>sys.stderr, "dagG.shape", dadG.shape
            #print >>sys.stderr, "momG.shape", momG.shape
            pedOut.write('\t'.join(dad + list(dadG))+'\n')
            pedOut.write('\t'.join(mom + list(momG))+'\n')
    
    for i,ch in enumerate(children):
        #print >>sys.stderr, "childG.shape", childG.shape
        pedOut.write('\t'.join(persons[ch] + list(childG[i]))+'\n')
    start = False

n = 0
cur_fam = None
family=[]
with open(htFn, 'r') as f:
    HD = {v:i for i,v in 
          enumerate(f.readline().strip('\n\r').split('\t'))}
    for l in f:
        fId, pId, htId, ht = l.strip('\n\r').split('\t')
        if cur_fam and fId != cur_fam:
            #print 'A'
            # process cur_fam
            processFamily(family)
            family=[]
        #print 'B'
        cur_fam = fId
        family.append({'fId':fId, 'pId':pId, 
                       'htId':htId, 'ht':ht})
        n += 1
        #if n == 80:
        #    break

# process the last family
processFamily(family);

bimOut = open(bimFn, 'w')
with open(posFn, 'r') as f:
    f.readline()
    for i,l in enumerate(f):
        cs = l.strip('\n\r').split('\t')
        cnt=defaultdict(int)

        for k,v in ALS[i].items():
            cnt[k] += v

        assert (('0' in cnt and len(cnt) < 4) 
                or (not '0' in cnt and len(cnt) < 3))

        if '0' in cnt and len(cnt) == 1:
            majorA = '0'
            minorA = '0'
        else:
            if '0' in cnt:
                del cnt['0']
            als = sorted(cnt.items(), key=lambda x: x[1])
        
            if len(als) == 1:
                minorA = '0'
                majorA = als[0][0]
            else:
                minorA = als[0][0]
                majorA = als[1][0]

        bimOut.write('\t'.join([cs[0], OBIM[(cs[0],cs[1])][1], '0', 
                                cs[1], minorA, majorA]) +'\n')

pedOut.close()
bimOut.close()

