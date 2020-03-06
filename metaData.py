#!/bin/env python

import os
#from pylab import *
import sys
import numpy as np

RLDR = {"mother":"mom", "father":"dad", "prb":"prb", "sibling":"sib"}

class MetaData:
    def buildFromNucFams(md):
        for fd in list(md.nucFams.values()):
            for pd in fd.memberInOrder:
                pd.familyId = fd.familyId
        md.persons = {pd.personId:pd for fd in list(md.nucFams.values()) for pd in fd.memberInOrder}

        md.trios = {}
        for fid,fd in list(md.nucFams.items()):
            
            for pd in fd.memberInOrder[2:]:
                t = Family()
                t.familyId = pd.personId 
                t.atts = dict(fd.atts)

                mom,dad = fd.memberInOrder[0:2]
                t.memberInOrder = [mom,dad,pd] 
                md.trios[t.familyId] = t

                pd.atts['momPersonId'] = mom.personId
                pd.atts['dadPersonId'] = dad.personId

class Family:
    pass

class Person:
    pass

def load_meta_data(MDFN=None):
    if not MDFN:
        MDFN = os.environ['PROJECT_DIR'] + "/metaData.txt"
        
    if 'SSC' in os.environ['PROJECT_DIR']:
        #print 'A'
        parsers = [_parse_SSC2, _parse_SSC, _parse_AGRE, _parse_AGRE2, _parse_PTCD]
    else:
        #print 'B'
        parsers = [_parse_AGRE2, _parse_SSC, _parse_SSC2, _parse_AGRE, _parse_PTCD]

    for parser in parsers:
        try:
            return parser(MDFN)
        except:
            pass

    return None

def _parse_PTCD(MDFN): 
    mD = np.genfromtxt(MDFN,delimiter='\t',dtype=None,names=True, case_sensitive=True)

    md = MetaData()
    md.nucFams = {}
   
    rlMap = {"M":"mom", "F":"dad", "A":"prb", "U":"sib"} 

    for R in mD:
        p = Person()

        p.atts = { x: R[x] for x in R.dtype.names }
        p.personId = R["Gleeson_ID_Copy"]
        p.smId = str(R["smId"])
        p.bamF = os.environ['PROJECT_DIR'] + "/bams/" + str(p.smId) + ".bam"
        p.role = rlMap[R['Patient_Name']]
        p.gender = R["Sex"]

        fmId = p.personId.split("-")[0]
        if fmId not in md.nucFams:
            f = Family()
            f.familyId = fmId
            f.memberInOrder = []
            f.atts = {} 
            md.nucFams[fmId] = f

        md.nucFams[fmId].memberInOrder.append(p)

    rlOrd = {rl:rli for rli,rl in enumerate("mom dad prb sib".split())}
    for fd in list(md.nucFams.values()):
       fd.memberInOrder = sorted(fd.memberInOrder,key=lambda p: rlOrd[p.role]) 

    md.buildFromNucFams()
    return md

def _parse_SSC(MDFN): 
    #print "_parse_SSC"
    mD = genfromtxt(MDFN,delimiter='\t',dtype=None,names=True, case_sensitive=True)
    allSamples = []
    allSamplesByFamily = {}

    ref = None

    md = MetaData()
    md.nucFams = {}
    
    for mdR in mD:
        if mdR['ready'] != 1:
            continue

        f = Family()
        f.familyId = str(mdR['familyId'])
        #f.atts = { x: mdR[x] for x in mdR.dtype.names }
        f.atts = {}
        
        def toPerson(rl):
            p = Person()
            if not ':' in mdR["BAM_ID_" + rl]:
                return
            p.personId,p.smId  = mdR["BAM_ID_" + rl].split(":") 
            p.atts = {}
            p.atts['bamF'] = mdR["BAM_FL_" + rl]
            p.atts['aff'] = 0
            if rl == 'prb':
                p.atts['aff'] = 1
            # TODO rethink!!!
            p.atts['batch'] = p.atts['bamF'].split("/")[3] 
            p.gender = mdR["GNDR_" + rl] 
            p.bamF = p.atts['bamF']
            p.role = rl
            return p

        f.memberInOrder = [toPerson(rl) for rl in "mom dad prb sib".split()]
        md.nucFams[f.familyId] = f
    md.buildFromNucFams()
    return md


def _parse_AGRE(MDFN): 
    #print "AGRE"
    class Individual:
        pass

    class NuclearFamily:
        pass

    class BamFile:
        pass

    bamFiles = {}
    individuals = {}
    nuclearFamilies = {}

    sex2GenderD = {"Female":"F", "F":"F", "Male":"M", "M":"M", '1':'M', '2':'F', 2:'F', 1:'M'}

    CT  = genfromtxt(MDFN,delimiter='\t',dtype=None,names=True, case_sensitive=True)
    for R in CT:
        ind = Individual()
        #ind.atts = { x: R[x] for x in CT.dtype.names }
        ind.atts = { x: R[x] for x in ['sampleId', 'mom', 'dad', 'category', 'bamF'] }
        ind.indId = R['indId']
        ind.smId = R['sampleId']
        ind.gender = sex2GenderD[R['sex']]
        ind.role = RLDR[R['relation']]
        ind.affectedStatus = R['aff']
        ind.category = R['category']

        ind.mother = None
        ind.father = None

        ind.nucFamily = None

        ind.bamF = None

        if R['bamF']:
            ind.bamF = BamFile()
            ind.bamF.ind = ind
            ind.bamF.fileName = R['bamF'] 
            ind.bamF.smId = R['sampleId']
            ind.smId = R['sampleId']
            assert ind.bamF.smId not in bamFiles
            bamFiles[ind.bamF.smId] = ind.bamF

        assert ind.indId not in individuals
        individuals[ind.indId] = ind

    for ind in list(individuals.values()):
        if ind.atts['mom'] != "0" and ind.atts['mom'] in individuals:
            ind.mother = individuals[ind.atts['mom']]
        if ind.atts['dad'] != "0" and ind.atts['dad'] in individuals:
            ind.father = individuals[ind.atts['dad']]



    for ind in list(individuals.values()):
        if not (ind.father and ind.mother):
            continue
        ncFId = ind.mother.indId + "_" + ind.father.indId
        if ncFId in nuclearFamilies:
            ncF = nuclearFamilies[ncFId]
        else: 
            ncF = NuclearFamily()
            ncF.ncFId = ncFId
            ncF.mother = ind.mother
            ncF.father = ind.father
            ncF.children = []
            nuclearFamilies[ncFId] = ncF

        ncF.children.append(ind)

        ind.nucFamily = ncF 

    md = MetaData()
    md.nucFams = {}

    for ncFId,ncF in sorted(nuclearFamilies.items()):
        if not (ncF.mother.bamF and ncF.father.bamF):
            continue
        chld = sorted([ch for ch in ncF.children if ch.bamF],key=lambda x: x.indId)
        if not chld:
            continue

        mom = ncF.mother
        dad = ncF.father

        mmbrs = [mom, dad] + chld
        rls = ["mom", "dad"] + len(chld) * [None]

        f = Family()

        f.familyId = ncFId
        f.atts = {"category":",".join({ind.atts['category'] for ind in mmbrs})}

        def toPerson(ind,rl):
            p = Person()
            p.personId = ind.indId
            p.smId = ind.smId
            p.atts = dict(ind.atts)
            # TODO rethink!!!
            p.atts['batch'] = p.atts['bamF'].split("/")[4] 
            p.bamF = p.atts['bamF']
            p.gender = ind.gender
            if rl:
                p.role = rl
            else:
                p.role = "sib" if ind.affectedStatus==0 else "prb"
            return p

        f.memberInOrder = [toPerson(ind,rl) for ind,rl in zip(mmbrs,rls)]

        md.nucFams[f.familyId] = f

    md.buildFromNucFams()
    return md

def _parse_AGRE2(MDFN):
    #print 'AGRE2'

    class Individual:
        pass

    class NuclearFamily:
        pass

    class BamFile:
        pass

    bamFiles = {}
    individuals = {}
    nuclearFamilies = {}
    
    sex2GenderD = {"Female":"F", "F":"F", "Male":"M", "M":"M", '1':'M', '2':'F', 2:'F', 1:'M'}

    CT  = genfromtxt(MDFN,delimiter='\t',dtype=None,names=True, case_sensitive=True)

    for R in CT:
        ind = Individual()
        #ind.atts = { x: R[x] for x in CT.dtype.names }
        #ind.atts = { x: R[x] for x in ['sampleId', 'mom', 'dad', 'category', 'bamF'] }
        ind.atts = eval(R['atts'])
        ind.indId = R['indId']
        ind.smId = R['sampleId']
        ind.gender = sex2GenderD[R['gender']]
        #ind.role = RLDR[R['relation']]
        ind.role = R['role']
        #ind.affectedStatus = R['aff']
        if ind.role in ['mom', 'dad']:
            ind.affectedStatus = 0
        else:
            ind.affectedStatus = 0 if ind.role == 'sib' else 1
        #ind.category = R['category']
        ind.category = ind.atts['category']
        ind.mother = None
        ind.father = None

        ind.nucFamily = None

        ind.bamF = None

        if R['bamF']:
            ind.bamF = BamFile()
            ind.bamF.ind = ind
            ind.bamF.fileName = R['bamF'] 
            ind.bamF.smId = R['sampleId']
            ind.smId = R['sampleId']
            assert ind.bamF.smId not in bamFiles
            bamFiles[ind.bamF.smId] = ind.bamF

        assert ind.indId not in individuals
        individuals[ind.indId] = ind

    for ind in list(individuals.values()):
        if ind.atts['mom'] != "0" and ind.atts['mom'] in individuals:
            ind.mother = individuals[ind.atts['mom']]
        if ind.atts['dad'] != "0" and ind.atts['dad'] in individuals:
            ind.father = individuals[ind.atts['dad']]



    for ind in list(individuals.values()):
        if not (ind.father and ind.mother):
            continue
        ncFId = ind.mother.indId + "_" + ind.father.indId
        if ncFId in nuclearFamilies:
            ncF = nuclearFamilies[ncFId]
        else: 
            ncF = NuclearFamily()
            ncF.ncFId = ncFId
            ncF.mother = ind.mother
            ncF.father = ind.father
            ncF.children = []
            nuclearFamilies[ncFId] = ncF

        ncF.children.append(ind)

        ind.nucFamily = ncF 

    md = MetaData()
    md.nucFams = {}

    for ncFId,ncF in sorted(nuclearFamilies.items()):
        if not (ncF.mother.bamF and ncF.father.bamF):
            continue
        chld = sorted([ch for ch in ncF.children if ch.bamF],key=lambda x: x.indId)
        if not chld:
            continue
        
        mom = ncF.mother
        dad = ncF.father

        mmbrs = [mom, dad] + chld
        rls = ["mom", "dad"] + len(chld) * [None]

        f = Family()
        f.familyId = ncFId
        f.atts = {"category":",".join({ind.atts['category'] for ind in mmbrs})}

        def toPerson(ind,rl):
            p = Person()
            p.personId = ind.indId
            p.smId = ind.smId
            p.atts = dict(ind.atts)
            # TODO rethink!!!
            p.atts['batch'] = p.atts['bamF'].split("/")[4] 
            p.bamF = p.atts['bamF']
            p.gender = ind.gender
            if rl:
                p.role = rl
            else:
                p.role = "sib" if ind.affectedStatus==0 else "prb"
            p.atts['aff'] = 1 if p.role == 'prb' else 0
            return p

        f.memberInOrder = [toPerson(ind,rl) for ind,rl in zip(mmbrs,rls)]
        md.nucFams[f.familyId] = f

    md.buildFromNucFams()
    return md


def _parse_SSC2(MDFN):
    #print '_parse_SSC2'
    class Individual:
        pass

    class NuclearFamily:
        pass

    class BamFile:
        pass

    bamFiles = {}
    individuals = {}
    nuclearFamilies = {}

    sex2GenderD = {"Female":"F", "F":"F", "Male":"M", "M":"M", '1':'M', '2':'F', 2:'F', 1:'M'}

    CT  = genfromtxt(MDFN,delimiter='\t',dtype=None,names=True, case_sensitive=True)
    #print len(CT)
    for R in CT:
        #print '|' + R['bamF'] + '|'
        ind = Individual()
        #ind.atts = { x: R[x] for x in CT.dtype.names }
        #ind.atts = { x: R[x] for x in ['sampleId', 'mom', 'dad', 'category', 'bamF'] }
        ind.atts = eval(R['atts'])
        ind.indId = R['indId']
        ind.smId = R['sampleId']
        ind.gender = sex2GenderD[R['gender']]
        #ind.role = RLDR[R['relation']]
        ind.role = R['role']
        #ind.affectedStatus = R['aff']
        if ind.role in ['mom', 'dad']:
            ind.affectedStatus = 0
        else:
            ind.affectedStatus = 0 if ind.role == 'sib' else 1
        #ind.category = R['category']
        ind.category = ind.atts['category']
        ind.mother = None
        ind.father = None
        ind.nucFamily = None
        ind.bamF = None
        if R['bamF']:
            ind.bamF = BamFile()
            ind.bamF.ind = ind
            ind.bamF.fileName = R['bamF'] 
            ind.bamF.smId = R['sampleId']
            ind.smId = R['sampleId']
            assert ind.bamF.smId not in bamFiles
            bamFiles[ind.bamF.smId] = ind.bamF
        assert ind.indId not in individuals
        individuals[ind.indId] = ind
    #print len(individuals)
    for ind in list(individuals.values()):
        if ind.atts['mom'] != "0" and ind.atts['mom'] in individuals:
            ind.mother = individuals[ind.atts['mom']]
        if ind.atts['dad'] != "0" and ind.atts['dad'] in individuals:
            ind.father = individuals[ind.atts['dad']]



    for ind in list(individuals.values()):
        if not (ind.father and ind.mother):
            continue
        ncFId = ind.mother.indId.split('.')[0]
        if ncFId in nuclearFamilies:
            ncF = nuclearFamilies[ncFId]
        else: 
            ncF = NuclearFamily()
            ncF.ncFId = ncFId
            ncF.mother = ind.mother
            ncF.father = ind.father
            ncF.children = []
            nuclearFamilies[ncFId] = ncF

        ncF.children.append(ind)

        ind.nucFamily = ncF 

    md = MetaData()
    md.nucFams = {}

    for ncFId,ncF in sorted(nuclearFamilies.items()):
        if not (ncF.mother.bamF and ncF.father.bamF):
            continue
        chld = sorted([ch for ch in ncF.children if ch.bamF],key=lambda x: x.indId)
        if not chld:
            continue

        mom = ncF.mother
        dad = ncF.father

        mmbrs = [mom, dad] + chld
        rls = ["mom", "dad"] + len(chld) * [None]

        f = Family()
        f.familyId = ncFId
        f.atts = {"category":",".join({ind.atts['category'] for ind in mmbrs})}

        def toPerson(ind,rl):
            p = Person()
            p.personId = ind.indId
            p.smId = ind.smId
            p.atts = dict(ind.atts)
            # TODO rethink!!!
            #p.atts['batch'] = p.atts['bamF'].split("/")[4] 
            p.bamF = p.atts['bamF']
            p.gender = ind.gender
            if rl:
                p.role = rl
            else:
                p.role = "sib" if ind.affectedStatus==0 else "prb"
            p.atts['aff'] = 1 if p.role == 'prb' else 0
            return p

        f.memberInOrder = [toPerson(ind,rl) for ind,rl in zip(mmbrs,rls)]

        md.nucFams[f.familyId] = f

    md.buildFromNucFams()
    return md


def smIdPersonIdbamF(md):
    print('\t'.join('sampleId personId bamF'.split(' ')))
    for n in md.persons:
        print('\t'.join([ md.persons[n].smId, md.persons[n].personId, md.persons[n].bamF ] ))

def addZygosity(md):
    individuals = {}
    
    RMD="/gpfs/commons/home/iossifovi-488/AGRE_WG/metaDataRaw"
    catalog = RMD + "/batchC/150827-catalog.txt"
    small = RMD + "/batchB/NYGC.clean112915.txt"
    
    CF = open(catalog)
    cD = CF.readline().strip().split("\t")
    for ln,l in enumerate(CF):
        cs = l.strip("\n").split("\t",-1)
        assert( len(cs) == len(cD) )
        r = dict(list(zip(cD,cs)))
        individuals[r['Individual Code']] = r
    CF.close()
    
    CF = open(small)
    cD = CF.readline().strip().split("\t")
    for ln,l in enumerate(CF):
        cs = l.strip("\n").split("\t",-1)
        assert( len(cs) == len(cD) )
        r = dict(list(zip(cD,cs)))
        if r['indId'] not in individuals:
            if r['zygosity'] == '0':
                individuals[r['indId']] = {'Zygosity Type':' '}
            elif r['zygosity'] == '2':
                individuals[r['indId']] = {'Zygosity Type':'Twin'}
    CF.close()
    
    persons = md.persons
    
    for p,v in list(persons.items()):
        if p in individuals:
            v.atts['Zygosity Type'] = individuals[p]['Zygosity Type']
        else:
            v.atts['Zygosity Type'] = ''


