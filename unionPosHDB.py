#!/usr/bin/env python
from math import *
import numpy as np
import glob 
from time import time
import sys
from glob import glob

posFiles = []

if len(sys.argv) < 3:
    print "Usage: unionPosHDB.py HDB_1 .. HDB_n  UHDB"
    sys.exit(1)

if len(sys.argv)==3:
    # if there are only two argumantes the first is the directory name
    # with hdbs
    posFiles = sorted(glob(sys.argv[1] + '/*-pos.txt'))
else:
    for line in sys.argv[1:-1]:
        print line
        line = line.strip()
        posFiles.append(line + '-pos.txt')

UHDB = sys.argv[-1]

st = time()
dt = np.dtype([('Chromosome', 'S24'), ('Position', '<i8'), ('RefAllele', 'S1'), ('A', '<i4'), ('C', '<i4'), ('G', '<i4'), ('T', '<i4')])
POS = np.genfromtxt(posFiles[0],delimiter='\t',dtype=dt,names=True, case_sensitive=True)
print >>sys.stderr, "read POS", len(POS), "lines"
#cnts = np.array(POS[['A','C','G','T']].tolist())
cnts = np.vstack((POS['A'],POS['C'],POS['G'],POS['T']))
#cnts = array(POS[:][['A','C','G','T']])

for i,f in enumerate(posFiles[1:]):
    print >>sys.stderr, "processing file", i, time() -st
    POS1 = np.genfromtxt(f,delimiter='\t',dtype=dt,names=True, case_sensitive=True)
    print >>sys.stderr, "read POS1", len(POS1), "lines"

    assert(POS['Position'].all() == POS1['Position'].all())
    #cnts += np.array(POS1[['A','C','G','T']].tolist())
    cnts += np.vstack((POS1['A'],POS1['C'],POS1['G'],POS1['T']))
    #cnts += array(POS1[:][['A','C','G','T']])
    print >>sys.stderr, 'read family ', f, (time() - st) 

S = POS.shape
FPOS = open(UHDB + '-pos.txt', 'w') 
FPOS.write( '\t'.join('Chromosome Position RefAllele A C G T'.split(' ')) + '\n')
for n in range(S[0]):
    FPOS.write( '\t'.join(map(str, list(POS[n])[:3] + list(cnts[:,n]))) + '\n')
FPOS.close()
print >>sys.stderr, "Done"
