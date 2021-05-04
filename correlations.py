#!/usr/bin/env python

import os,sys,optparse
import numpy as np
from collections import Counter,defaultdict

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage=usage)

parser.add_option("-c", "--chrom", dest="chr", default="22",
                  metavar="CHROM", help="chromosome", type="string")
parser.add_option("-w", "--whites", dest="whites", default=0,
                  metavar="whites", help="only whites", type=int)
parser.add_option("-f", "--freq_min", dest="freq_min", default=0.001,
                  metavar="freq_min", help="min MAF included", type="float")
parser.add_option("-F", "--freq_max", dest="freq_max", default=1.0,
                  metavar="freq_max", help="max MAF excluded", type="float")
parser.add_option("-a", "--anti", dest="anti", default=0,
                  metavar="ANTI", help="anti child: 0|1", type="int")
parser.add_option("-r", "--result_dir", dest="res_dir", default="CORR",
                  metavar="RES_DIR", help="results directory name", type="string")
parser.add_option("-p", "--prefix", dest="prefix", default="CORR",
                  metavar="prefix", help="names of sets to be correlated", type="string")
parser.add_option("-S", "--famStruct", dest="famStructFn", default="quadsStruct.txt",
                  metavar="famStrucFn", help="family structure file", type="string")


ox, args = parser.parse_args()
print(ox, file=sys.stderr)

chr = ox.chr
freq_min = ox.freq_min
freq_max = ox.freq_max
anti = ox.anti
res_dir = ox.res_dir
W = ox.whites
pref = ox.prefix.split(',')
famStructFn = ox.famStructFn

# c is array of random 1,-1 with the number of rows
# equal to the number probands in quads
# and the number of columns equal to number of iterations (10000)
# it is precomputed and stored so that
# the same matrix is used for all chromosomes

gnpFn = pref[0] + '.npz'
gncFn = pref[1] +'.npz'
pos1Fn = pref[0]  + '.pos'
pos2Fn = pref[1]  + '.pos'
chFn = pref[1] + '.indv'
prFn = pref[0] + '.indv'
#famStructFn = 'goodFamStruct.txt'
#famStructFn = 'quadsStruct.txt'

famStruct = {}
with open(famStructFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        famStruct[cs[0]] = cs
    
pos1 = {}
with open(pos1Fn, 'r') as f:
    for i,l in enumerate(f):
        cs = l.strip('\n\r').split('\t')
        if 'chr' in cs[0]:
            cs[0] = cs[0][3:]
        pos1[tuple(map(int,cs[:2]))] = i

pos2 = {}
with open(pos2Fn, 'r') as f:
    for i,l in enumerate(f):
        cs = l.strip('\n\r').split('\t')
        if 'chr' in cs[0]:
            cs[0] = cs[0][3:]
        pos2[tuple(map(int,cs[:2]))] = i



# they are in mo,fa,... order 
parents = []
with open(prFn, 'r') as f:
    for l in f:
        parents.append( l.strip('\n\r') )

children = []
with open(chFn, 'r') as f:
    for l in f:
        if 'spark' in pref[1]:
            children.append( famStruct[l.strip('\n\r')][2] )
        else:
            children.append(l.strip('\n\r'))

f = np.load(gnpFn, 'rb')
print("Keys: %s" % list(f.keys()), file=sys.stderr)
gnp = f.get('GN').astype(float)
print("reading parents npz:", file=sys.stderr)
print("gnp.shape", gnp.shape, file=sys.stderr)
f.close()


f = np.load(gncFn, 'rb')
print("Keys: %s" % list(f.keys()), file=sys.stderr)
gnc = f.get('GN').astype(float)
print("reading children npz:", file=sys.stderr)
print("gnc.shape", gnc.shape, file=sys.stderr)
f.close()

prs = gnp
chn = gnc

corr = np.corrcoef(chn)
id = np.where(corr > 0.95)
l = len(id[0])
idd = set([id[0][i] for i in range(l) if id[0][i] != id[1][i]])

for i in idd:
    print(children[i])
