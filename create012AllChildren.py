#!/usr/bin/env python

import sys, os, gzip
import numpy as np
from collections import Counter,defaultdict


if len(sys.argv) < 5:
    print("Usage: create012AllChildren.py <HDB> <out-012-file> phased=1 or 0 <famNum>")
    exit(1)

hdb = sys.argv[1]
out = sys.argv[2]
phased = bool(sys.argv[3])
famNum = int(sys.argv[4])

print(hdb, out, phased, file=sys.stderr)

badGens = set(list('oucd'))
ambiguous = {
    'AC':'M', 
    'AG':'R', 
    'AT':'W', 
    'CG':'S', 
    'CT':'Y', 
    'GT':'K'
}
ramb = {v:k for k,v in list(ambiguous.items())}
for t,c in list(ambiguous.items()):
    ambiguous[t[1] + t[0]] = c

class P:
    def __init__(p,cs):
        p.ch,p.pos,p.refA,p.A,p.C,p.G,p.T = cs
        p.pos,p.A,p.C,p.G,p.T = list(map(int,[p.pos,p.A,p.C,p.G,p.T]))
        aad = {a:n for a,n in zip("ACGT",[p.A,p.C,p.G,p.T])} 
        ad = {a:n for a,n in list(aad.items()) if a!=p.refA and n>0}
        p.AI = {p.refA:0}
        for ai,(a,n) in enumerate(sorted(list(ad.items()),key=lambda x: (-x[1],x[0]))):
            p.AI[a] = ai+1

        p.unphasedD = {}
        p.phasedD = {}
        for a,ai in list(p.AI.items()):
            p.unphasedD[a + a] = "%d/%d" % (ai,ai)
            p.phasedD[a + a] = "%d|%d" % (ai,ai)
            for b,bi in list(p.AI.items()):
                if a==b: continue
                fi,si = sorted([ai,bi])
                
                p.unphasedD[a + b] = "%d/%d" % (fi,si)
                p.phasedD[a + b] = "%d|%d" % (ai,bi)
                p.unphasedD[ambiguous[a + b]*2] = "%d/%d" % (fi,si)
                p.phasedD[ambiguous[a + b]*2] = "%d/%d" % (fi,si)

        # unphasedD {"AA": "0/0", "AC": 0/1, "CA": 0/1, "RR": "0/1"}
        # phasedD   {"AA": "0|0", "AC": 0|1, "CA": 1|0, "RR": "0/1"}
        
       

    def vcfGen(p,AT,ANT, phased):
        if AT in badGens:
            return '-1'
        D = p.phasedD if phased else p.unphasedD
        #return D[AT+ANT]
        return sum(map(int,[D[AT+ANT][0], D[AT+ANT][2]]))
    
        # if phased==True:  0|0  1|0 0|1, ...
        # else 0/0, 0/1, 0/2
        # * if AT or ANT are small letters: ./. or .|. 
        # * if AT is a two-nucleotide code, assert that AT == ANT, and 
        # output unphased regardless of phased parameter

POS = []
posFn = hdb + '-pos.txt'
with open(posFn, 'r') as f:
    f.readline()

    for l in f:
        POS.append(P(l.strip('\n\r').split('\t')))

FM_CH = {} 
FM_HPTH = defaultdict(dict) 

hpthFn = hdb + '-hpth.txt'

if os.path.isfile(hpthFn):
    f = open(hpthFn)
elif  os.path.isfile(hpthFn + '.gz'):
    f = gzip.open(hpthFn + '.gz')
else:
    print("no file ", hpthFn + '[.gz]')
    exit(1)

f.readline()
for li,l in enumerate(f):
    if li % 100 == 0:
        print(li, "...")
    fmId,perId,thrT,thr = l.strip('\n\r').split('\t')
    if thrT == 'MT' or thrT == 'FT':
            FM_HPTH[perId][thrT] = thr
    #if li == 39:
    #    break
f.close()

NH ,= set(len(x) for x in list(FM_HPTH.values()) )
assert NH == 2

perIdOrder = sorted(FM_HPTH.keys())

### TODO: get that somehow!!!!
"""
with open(faiFn, 'r') as f:
    chrLen = {'chr'+str(x):0 for x in range(1,23)}
    for l in f:
        cs = l.strip('\n\r').split('\t')
        if cs[0] in chrLen:
            chrLen[cs[0]] = int(cs[1])
"""        
OF = open(out,"w")
# OF = sys.stdout 

#OF.write('##fileformat=VCFv4.2\n')
#OF.write('##fileDate=20171012\n')
#OF.write('##source=createVcfAllParents.py\n')
#for ch in sorted({p.ch for p in POS}):
#    OF.write('##contig=<ID=%s,length=%d>\n' % (ch,chrLen[ch]))
#OF.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')


OF.write("\t".join('#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT'.split(' ') + 
            [perId for perId in perIdOrder]) + "\n")

filterStats = defaultdict(int)
for pi,p in enumerate(POS):
    #if (pi % 1000) == 0:
    #    print "writing out position", pi, "..."
    """
    # no laternative alleles
    if len(p.AI) == 1:
        filterStats['no_alt'] += 1
        print >>sys.stderr, 'no_alt', p.ch, p.pos
        continue
    
    # more than one alternative 
    if len(p.AI) > 2:
        filterStats['more_alt'] += 1
        print >>sys.stderr, 'more_alt', p.ch, p.pos
        continue
    
    # refA not seen
    cntsA = dict(zip('ACGT',[p.A,p.C,p.G,p.T]));
    if cntsA[p.refA] == 0:
        filterStats['no_ref'] += 1
        print >>sys.stderr, 'no_ref', p.ch, p.pos
        continue
    
    # less than 90% genotyped
    # this is hack, assumes that there two children in each family
    NG = sum(cntsA.values())
    if NG < 0.9 * famNum * 2 * 2:
        filterStats['no_ref'] += 1
        print >>sys.stderr, 'less_then_90%', p.ch, p.pos, NG, famNum
        continue
    """
    #if 100.*(NG - cntsA[p.refA])/NG < th:
    #    continue

    cs = [p.ch,
          str(p.pos), 
          "%s:%d" % (p.ch,p.pos), 
          p.refA, 
          ",".join([a for a,ai in sorted(list(p.AI.items()),key=lambda x: -x[1]) if ai>0]), 
          ".", 
          ".", 
          ".", 
          "GT"]
    for fmId in perIdOrder:
        fmD = FM_HPTH[fmId]
        cs += [p.vcfGen(fmD['MT'][pi],fmD['FT'][pi],phased)] 
    OF.write("\t".join(map(str,cs)) + "\n") 
    # if pi > 100000:
    #     break
OF.close()
