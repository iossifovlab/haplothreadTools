#!/usr/bin/env python

import os,sys
import numpy as np


if len(sys.argv) <2:
    print("Usage: reformat012.py <prefix: parents|children>")
    exit(1)
    

pref = sys.argv[1]
fin = pref + '.012'
fout = pref + '.npz'
dname = pref

indFn = pref + '.indv'
posFn = pref + '.pos'

N = 1100000

pos = []
INDV = open(indFn, 'w')
with open(fin, 'r') as f:
    cs = f.readline().strip('\n\r').split('\t')
    for p in cs[9:]:
        INDV.write(p+'\n')
    INDV.close()
    gn = np.empty((N, len(cs[9:])), dtype='int')
    gn[:] = 3
    for k,l in enumerate(f):
        cs = l.strip('\n\r').split('\t')
        gn[k,:] = cs[9:]
        chrom, p = cs[:2]
        #pos.append(['chr'+chrom, p])
        pos.append([chrom, p])
        #if k == 100:
        #    break

id = (gn[:,0] != 3)
GN = gn[id,:].transpose()

d_file = open(fout, 'wb')
np.savez(d_file, GN=GN)
d_file.close()

with open(posFn, 'w') as f:
    for p in pos:
        f.write('\t'.join(p) +'\n')


f = np.load(fout, 'rb')
print("Keys: %s" % list(f.keys()), file=sys.stderr)
gnp = f.get('GN')
print("reading parents npz:", file=sys.stderr)

# this can be commented

for n in range(gnp.shape[0]):
    a = ''.join(GN[n,:].astype(str))
    b = ''.join(gnp[n,:].astype(str))
    assert(a==b)

print("Finished verification", file=sys.stderr)


