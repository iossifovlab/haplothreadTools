#!/usr/bin/env python

import sys, os, argparse
from collections import Counter
import numpy as np
from metaData import *
import gzip

# modes for filtering by position:
#    1. chr:begChrPos-endChrPos
#    2. begPosIndex endPosIndex  (TODAY!!!)
#    3. list of specific positions:
#        3.1. comma separated ch1:pos1,ch2:pos2,...
#        3.2. bed file for specific positions
#    4. frequency based:
#        all common (altF larger that a given parameters)
#        all rare (altF smaller that a given parameters)
#        all ultra-rare (altF smaller that a given parameters)
#        all non-ultra-rare (altF smaller that a given parameters)
#
# modes for filtering by haplothread
#   1. by familyId(s) (from comma-separated parameter or from a file)
#   2. by childId(s) (from comma-separated parameter or from a file)
#   3. by phenotypic property??? 
#   4. by haplothread type (MT,MNT,FT,FNT) 
#   NOTE: that the haplothread subset should update the frequency in the position file.


'''
if len(sys.argv) < 5:
    print "Usage: posSubsetHDB.py IHDB begPosIndex endPosIndex OHDB"
    sys.exit(1)
'''
md = None
ht2cnt = {
    'M':[1, 1, 0, 0],
    'R':[1, 0, 1, 0],
    'W':[1, 0, 0, 1],
    'S':[0, 1, 1, 0],
    'Y':[0, 1, 0, 1],
    'K':[0, 0, 1, 1],
    'A':[2, 0, 0, 0],
    'C':[0, 2, 0, 0],
    'G':[0, 0, 2, 0],
    'T':[0, 0, 0, 2],
    'c':[0, 0, 0, 0],
    'd':[0, 0, 0, 0],
    'o':[0, 0, 0, 0],
    'u':[0, 0, 0, 0]
}

def posIndexSubset(IHDB, OHDB, B,E):
    B = int(B)
    E = int(E)

    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())
    n = -1 
    for l in IP:
        n+=1
        if n<B:
            continue
        if n>=E:
            break
        OP.write(l)
    IP.close()
    OP.close()

    OH = open(OHDB+"-hpth.txt","w")
    IH = gzip.open(IHDB+"-hpth.txt.gz")

    OH.write(IH.readline())

    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4

        cs[3] = cs[3][B:E]
        OH.write("\t".join(cs)+"\n")
    IH.close()

def posRegionSubset(IHDB, OHDB, ch, B, E):
    B = int(B)
    E = int(E)

    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())
    n = -1 
    BB = False
    EE = False
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        if not cs[0] == ch:
            print("Wrong chromosome for IHDB: " + cs[0] + "!=" + ch) 
            exit(1)
        n+=1
        if int(cs[1])<B:
            continue
        if int(cs[1])>=E:
            if not EE:
                EE = n-1
            break
        if not BB:
            BB = n
        OP.write(l)
    IP.close()
    OP.close()

    B = BB
    if EE:
        E = EE
    OH = open(OHDB+"-hpth.txt","w")
    IH = gzip.open(IHDB+"-hpth.txt.gz")

    OH.write(IH.readline())

    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        if EE == False:
            cs[3] = cs[3][B:]
        else:
            cs[3] = cs[3][B:E]
        OH.write("\t".join(cs)+"\n")
    IH.close()

def posPositionsSubset(IHDB, OHDB, pos):
    # the commentd code assumed sorted pos,
    # this may not be the case with mapping to the new reference
    # so I change this code
    print("len(pos): ",len(pos), file=sys.stderr)
    """
    indx = []
    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())
    L = len(pos)
    k = 0
    n = -1 
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        n+=1
        if int(cs[1])<pos[k]:
            continue
        elif int(cs[1]) == pos[k]:
            OP.write(l)
            indx.append(n)
            k+=1
        if k == L:
            break
    IP.close()
    OP.close()
    print >>sys.stderr, "len(indx): ",len(indx)
    """
    indx = []
    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())
    L = len(pos)
    POS = [l.strip("\n\r").split("\t") for l in IP]
    ipos = {v[1]:[i]+v for i,v in enumerate(POS)}
    for p in pos:
        OP.write('\t'.join(ipos[str(p)][1:])+'\n')
        indx.append(ipos[str(p)][0])
    IP.close()
    OP.close()
    print("len(indx): ",len(indx), file=sys.stderr)
        
    OH = open(OHDB+"-hpth.txt","w")
    IH = gzip.open(IHDB+"-hpth.txt.gz")

    OH.write(IH.readline())

    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        cs[3] = "".join([cs[3][k] for k in indx])
        OH.write("\t".join(cs)+"\n")
    IH.close()


def posMapSubset(IHDB, OHDB, maps):
    # mapping from hg19 to hg38
    # maps has hg38 in the first column (sorted)
    # hg19 is in the second column

    print("len(maps): ",len(maps), file=sys.stderr)
    indx = []
    # the key will be hg19
    MAPS = {v[1]:v[0] for v in maps}
    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())
    L = len(MAPS)
    POS = [l.strip("\n\r").split("\t") for l in IP]
    ipos = {v[1]:[i]+v for i,v in enumerate(POS)}
    for p in sorted(MAPS, key=lambda x: int(MAPS[x])):
        OP.write('\t'.join([ipos[p][1]]+[MAPS[p]]+ipos[p][3:])+'\n')
        indx.append(ipos[p][0])
    IP.close()
    OP.close()
    print("len(indx): ",len(indx), file=sys.stderr)
        
    OH = open(OHDB+"-hpth.txt","w")
    IH = gzip.open(IHDB+"-hpth.txt.gz")

    OH.write(IH.readline())

    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        cs[3] = "".join([cs[3][k] for k in indx])
        OH.write("\t".join(cs)+"\n")
    IH.close()

def posBedFileSubset(IHDB, OHDB, bedFile):
    BF = open(bedFile, "r")
    bed = []
    for l in BF:
        bed.append(l.strip("\n\r").split("\t"))
    L = len(bed)
    indx = []
    k = 0
    B = int(bed[k][1])
    E = int(bed[k][2])

    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())
    n = -1 
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        if not cs[0] == bed[k][0]:
            print("Wrong chromosome for IHDB: " + cs[0] + "!=" + ch) 
            exit(1)
        n+=1
        if int(cs[1])<B:
            continue
        elif int(cs[1]) <= E:
            OP.write(l)
            indx.append(n)
        else:
            k+=1
            if k == L:
                break
            B = int(bed[k][1])
            E = int(bed[k][2])
            
    IP.close()
    OP.close()

    OH = open(OHDB+"-hpth.txt","w")
    IH = gzip.open(IHDB+"-hpth.txt.gz")

    OH.write(IH.readline())

    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        cs[3] = "".join([cs[3][k] for k in indx])
        OH.write("\t".join(cs)+"\n")
    IH.close()

def posChrPosFileSubset(IHDB, OHDB, chrPosFile):
    BF = open(chrPosFile, "r")
    bed = {}

    for l in BF:
        cs = l.strip("\n\r").split("\t") 
        bed[(cs[0],cs[1])] = 1
    L = len(bed)
    indx = []
    #k = 0
    #B = int(bed[k][1])
    #E = int(bed[k][2])

    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())
    n = -1 
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        n+=1
        if not (cs[0],cs[1]) in bed:
            continue
        OP.write(l)
        indx.append(n)
    IP.close()
    OP.close()

    OH = open(OHDB+"-hpth.txt","w")
    if os.path.isfile(IHDB+"-hpth.txt.gz"):
        IH = gzip.open(IHDB+"-hpth.txt.gz")
    elif os.path.isfile(IHDB+"-hpth.txt"):
        IH = open(IHDB+"-hpth.txt")
    else:
        print("Cannot open file", IHDB+"-hpth.txt*", file=sys.stderr)

    l = IH.readline()
    type_l = type(l)
    if type_l == bytes:
        OH.write(l.decode('utf8'))
    else:
        OH.write(l)

    for l in IH:
        if type_l:
            cs = l.decode('utf8').strip("\n\r").split("\t")
        else:
            cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        cs[3] = "".join([cs[3][k] for k in indx])
        OH.write("\t".join(cs)+"\n")
    IH.close()

def posMinAltAFreqSubset(IHDB, OHDB, minAltAFreq):
    MAP = {v:i for i,v in enumerate(list('ACGT'))}
    indx = []
    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())   
    n = -1 
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        n+=1
        cnts = list(map(float, cs[3:]))
        total = sum(cnts)
        altAFreq = 0
        if total > 0:
            altAFreq = (total - cnts[MAP[cs[2]]])/total
        if altAFreq < minAltAFreq:
            continue
        OP.write(l)
        indx.append(n)
    IP.close()
    OP.close()

    OH = open(OHDB+"-hpth.txt","w")
    IH = gzip.open(IHDB+"-hpth.txt")

    OH.write(IH.readline())

    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        cs[3] = "".join([cs[3][k] for k in indx])
        OH.write("\t".join(cs)+"\n")
    IH.close()

def posMaxAltAFreqSubset(IHDB, OHDB, maxAltAFreq):
    MAP = {v:i for i,v in enumerate(list('ACGT'))}
    indx = []
    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())   
    n = -1 
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        n+=1
        cnts = list(map(float, cs[3:]))
        total = sum(cnts)
        altAFreq = 0
        if total > 0:
            altAFreq = (total - cnts[MAP[cs[2]]])/total
        if altAFreq > maxAltAFreq:
            continue
        OP.write(l)
        indx.append(n)
    IP.close()
    OP.close()

    OH = open(OHDB+"-hpth.txt","w")
    IH = gzip.open(IHDB+"-hpth.txt")

    OH.write(IH.readline())

    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        cs[3] = "".join([cs[3][k] for k in indx])
        OH.write("\t".join(cs)+"\n")
    IH.close()

def posUltraRareSubset(IHDB, OHDB):
    MAP = {v:i for i,v in enumerate(list('ACGT'))}
    indx = []
    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())   
    n = -1 
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        n+=1
        cnts = list(map(float, cs[3:]))
        total = sum(cnts)
        altAFreq = 0
        if total > 0:
            altAFreq = (total - cnts[MAP[cs[2]]])
        if not altAFreq == 1:
            continue
        OP.write(l)
        indx.append(n)
    IP.close()
    OP.close()

    OH = open(OHDB+"-hpth.txt","w")
    IH = gzip.open(IHDB+"-hpth.txt")

    OH.write(IH.readline())

    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        cs[3] = "".join([cs[3][k] for k in indx])
        OH.write("\t".join(cs)+"\n")
    IH.close()

def posNucFamIdsSubset(IHDB, OHDB, famIds):
    OH = open(OHDB+"-hpth.txt","w")
    if os.path.isfile(IHDB+"-hpth.txt.gz"):
        IH = gzip.open(IHDB+"-hpth.txt.gz")
    elif os.path.isfile(IHDB+"-hpth.txt"):
        IH = open(IHDB+"-hpth.txt")
    else:
        print("Cannot open file", IHDB+"-hpth.txt*", file=sys.stderr)

    OH.write(IH.readline())
    st = True
    prev = None
    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        if not cs[0] in famIds:
            continue
        # to create pos file we have to allocate space based in the length of haplo-theread
        # do it only when encountered the first family to be included
        # then increment counts only for the first child in the family (four lines)
        if st:
            counts = np.array([[0,0,0,0]]*len(cs[3]))
            st = False
        if prev == None or not prev == cs[0]:
            prev = cs[0]
            k = -1
        k +=1
        if k < 4:
            counts += [ht2cnt[i] for i in cs[3]]
        OH.write(l)
    IH.close()
    OH.close()

    counts /= 2 # we double counts in ht2cnt to make all counts integer

    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())   
    n = -1 
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        n+=1
        OP.write('\t'.join(cs[:3] + list(map(str, counts[n,:]))) + '\n')
    IP.close()
    OP.close()

def posChildIdsSubset(IHDB, OHDB, childIds):
    OH = open(OHDB+"-hpth.txt","w")
    if os.path.isfile(IHDB+"-hpth.txt.gz"):
        IH = gzip.open(IHDB+"-hpth.txt.gz")
    elif os.path.isfile(IHDB+"-hpth.txt"):
        IH = open(IHDB+"-hpth.txt")
    else:
        print("Cannot open file", IHDB+"-hpth.txt*", file=sys.stderr)

    OH.write(IH.readline())
    st = True
    prev = None
    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        if not cs[1] in childIds:
            continue
        # to create pos file we have to allocate space based in the length of haplo-theread
        # do it only when encountered the first family to be included
        # then increment counts only for the first child in the family (four lines)
        if st:
            counts = np.array([[0,0,0,0]]*len(cs[3]))
            st = False
        if prev == None or not prev == cs[1]:
            prev = cs[1]
            k = -1
        k +=1
        if k < 4:
            counts += [ht2cnt[i] for i in cs[3]]
        OH.write(l)
    IH.close()
    OH.close()

    counts /= 2 # we double counts in ht2cnt to make all counts integer

    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())   
    n = -1 
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        n+=1
        OP.write('\t'.join(cs[:3] + list(map(str, counts[n,:]))) + '\n')
    IP.close()
    OP.close()

def posHaploThreadSubset(IHDB, OHDB, haplothread):
    OH = open(OHDB+"-hpth.txt","w")
    IH = gzip.open(IHDB+"-hpth.txt")

    OH.write(IH.readline())
    st = True
    prev = None
    for l in IH:
        cs = l.strip("\n\r").split("\t")
        assert len(cs)==4
        if not cs[2] == haplothread:
            continue
        # to create pos file we have to allocate space based in the length of haplo-theread
        # do it only when encountered the first family to be included
        # then increment counts only for the first child in the family (four lines)
        if st:
            counts = np.array([[0,0,0,0]]*len(cs[3]))
            st = False
        if prev == None or not prev == cs[0]:
            prev = cs[0]
            k = -1
        k +=1
        if k < 4:
            counts += [ht2cnt[i] for i in cs[3]]  
            # this should be reworked if parent transmits/not transmitts different alleles to different children
            # by now we keep only what is transmitted/not transmitted to the first child in the family
        OH.write(l)
    IH.close()
    OH.close()

    counts /= 2 # we double counts in ht2cnt to make all counts integer

    OP = open(OHDB+"-pos.txt","w")
    IP = open(IHDB+"-pos.txt");
    OP.write(IP.readline())   
    n = -1 
    for l in IP:
        cs = l.strip("\n\r").split("\t")
        n+=1
        OP.write('\t'.join(cs[:3] + list(map(str, counts[n,:]))) + '\n')
    IP.close()
    OP.close()


def file_len(fname):
    with open(fname) as f:
        for i,l in enumerate(f):
            pass
    return i

def main():

   #usage="usage: %prog [options]"
   
    parser = argparse.ArgumentParser(prog='posSubsetHDB.py')

    parser.add_argument("IHDB",  \
        default=None, action='store',
        metavar="IHDB", help="IHDB: input DB" )
    parser.add_argument("OHDB", \
            default=None, action='store',
            metavar="OHDB", help="OHDB: inoutput DB" )
    parser.add_argument("-r", "--region", dest="region", \
	     default=None, action='store',
             metavar="region", help="region chr[:begChrPos-endChrPos]" )
    parser.add_argument("-i", "--indexPos", dest="indexPos", \
	     default=None, action='store',
             metavar="indexPos", help="beginIndex:endIndex" )
    parser.add_argument("-p", "--positions", dest="positions", \
	     default=None, action='store',
             metavar="positions", help="comma separated pos1,pos2,.. or fileName" ) # ignore ch1:pos1,ch2:pos2,...
    parser.add_argument("-m", "--map_positions", dest="maps", \
	     default=None, action='store',
             metavar="maps", help="file containing newRefPos and oldRefPos" )
    parser.add_argument("-b", "--bedFile", dest="bedFile", \
	     default=None, action='store',
             metavar="bedFile", help="bed file name" )
    parser.add_argument("-B", "--chrPosFile", dest="chrPosFile", \
	     default=None, action='store',
             metavar="chrPosFile", help="chrom pos file name" )    
    parser.add_argument("-c", "--common", dest="minAltAFreq", \
	     default=None, action='store',
             metavar="minAltAFreq", help="altA frequency larger that a given parameter" )
    parser.add_argument("-F", "--maxAltAFreq", dest="maxAltAFreq", \
	     default=None, action='store',
             metavar="maxAltAFreq", help="altA frequency less that a given parameter" )
    parser.add_argument("-u", "--ultraRare", dest="ultraRare", \
             default=False, action='store_true', 
             help="only one alternative case in all collection" )
    parser.add_argument("-n", "--nucFamIds", dest="nucFamIds", \
	     default=None, action='store',
             metavar="nucFamIds", help="from comma-separated parameter or from a file" )
    parser.add_argument("-k", "--childIds", dest="childIds", \
	     default=None, action='store',
             metavar="childIds", help="from comma-separated parameter or from a file" )
    parser.add_argument("-P", "--phenotype", dest="phenotype", \
	     default=None, action='store',
             metavar="phenotype", help="project_dir:prb or project_dir:sib; more will be defined later" )
    parser.add_argument("-T", "--haploThread", dest="haploThread", \
	     default=None, action='store',
             metavar="haploThread", help="comma separated sublist of MT,MNT,FT,FNT" )

    args = parser.parse_args()
    print(args.IHDB)
    print(args.OHDB)
    #print args.indexPos

    #if args.region == None and args.indexPos == None:
    #    print "Number of positions: ", file_len(args.IHDB + '-pos.txt')
    #    print "Number of haplo-threads: ", file_len(args.IHDB + '-hpth.txt')

    if args.indexPos:
        B,E = args.indexPos.split(':')
        posIndexSubset(args.IHDB,args.OHDB,B,E)

    if args.region:
        ch, p = args.region.split(':')
        if "-" in p:
            B,E = p.split('-')
        else:
            B = p
            E = 300000000
        posRegionSubset(args.IHDB, args.OHDB, ch, B, E)
 
    if args.positions:
        if "," in args.positions:
            positions = sorted(map(int,args.positions.split(",")))
        elif os.path.isfile(args.positions):
            F = open(args.positions,"r")
            positions = []
            for l in F:
                positions.append(int(l.strip("\n\r")))
            F.close()
        else:
            print("file does not exist: |" + args.positions + "|")
            exit(1)
        posPositionsSubset(args.IHDB, args.OHDB, positions)

    if args.maps:
        if os.path.isfile(args.maps):
            with open(args.maps,"r") as F:
                maps = [l.strip("\n\r").split("\t") for l in F]
        else:
            print("file does not exist: |" + args.positions + "|")
            exit(1)
        posMapSubset(args.IHDB, args.OHDB, maps)
        
    if args.bedFile:
        posBedFileSubset(args.IHDB, args.OHDB, args.bedFile)

    if args.chrPosFile:
        posChrPosFileSubset(args.IHDB, args.OHDB, args.chrPosFile)        

    if args.minAltAFreq:
        posMinAltAFreqSubset(args.IHDB, args.OHDB, float(args.minAltAFreq))
        
    if args.maxAltAFreq:
        posMaxAltAFreqSubset(args.IHDB, args.OHDB, float(args.maxAltAFreq))

    if args.ultraRare:
        posUltraRareSubset(args.IHDB, args.OHDB)

    if args.nucFamIds:
        if "," in args.nucFamIds:
            famIds = args.nucFamIds.split(",")
        elif os.path.isfile(args.nucFamIds):
            F = open(args.nucFamIds,"r")
            famIds = []
            for l in F:
                famIds.append(l.strip("\n\r"))
            F.close()
        else:
            print("file does not exist: |" + args.nucFamIds + "|")
            exit(1)
        posNucFamIdsSubset(args.IHDB, args.OHDB, famIds)
                
    if args.childIds:
        if "," in args.childIds:
            childIds = args.childIds.split(",")
            if len(childIds[-1]) == 0:
                childIds = childIds[:-1]
        elif os.path.isfile(args.childIds):
            F = open(args.childIds,"r")
            childIds = []
            for l in F:
                childIds.append(l.strip("\n\r"))
            F.close()
        else:
            print("file does not exist: |" + args.childIds + "|")
            exit(1)
        posChildIdsSubset(args.IHDB, args.OHDB, childIds)

    if args.phenotype:
        if not ":" in args.phenotype:
            print("malformed phenotype: should be project_dir:prb or project_dir:sib")
            exit(1)
        project, phenotype = args.phenotype.split(":")
        if not phenotype in ["prb", "sib"]:
            print("invalid phenotype tipe, should be 'prb' or 'sib': " + "|" + phenotype + "|")
            exit(1)
        md = load_meta_data(project + "/metaData.txt")
        childIds =  [member.personId for fam in list(md.nucFams.values()) for member in fam.memberInOrder[2:] if member.role == phenotype]
        posChildIdsSubset(args.IHDB, args.OHDB, childIds)

    if args.haploThread:
        if not args.haploThread in "MT,MNT,FT,FNT":
            print("haploThread should be on of 'MT,MNT,FT,FNT'" + "|" + args.haploThread + "|")
            exit(1)
        posHaploThreadSubset(args.IHDB, args.OHDB, args.haploThread)                

if __name__ == "__main__":
	main()
