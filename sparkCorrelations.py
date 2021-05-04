#!/usr/bin/env python

import os,sys,optparse
import numpy as np
from collections import Counter,defaultdict

usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage=usage)

parser.add_option("-p", "--prefix", dest="prefix", default="CORR",
                  metavar="prefix", help="names of sets to be correlated", type="string")

ox, args = parser.parse_args()
print(ox, file=sys.stderr)

pref = ox.prefix.split(',')

# c is array of random 1,-1 with the number of rows
# equal to the number probands in quads
# and the number of columns equal to number of iterations (10000)
# it is precomputed and stored so that
# the same matrix is used for all chromosomes

gnpFn = pref[0] + '.npz'
gncFn = pref[1] +'.npz'
chFn = pref[1] + '.indv'
prFn = pref[0] + '.indv'

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
