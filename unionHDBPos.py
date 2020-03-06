#!/usr/bin/env python
import numpy as np
from time import time
import sys,os
import h5py
from glob import glob

posFiles = []
hpthFiles = []

if len(sys.argv) < 3:
    print("Usage: unionHDB.py HDB_1 .. HDB_n  UHDB")
    sys.exit(1)

if len(sys.argv)==3:
    # if there are only two argumantes the first is the directory name
    # with hdbs
    posFiles = sorted(glob(sys.argv[1] + '/*-pos.h5'))
else:
    for line in sys.argv[1:-1]:
        print(line, file=sys.stderr)
        line = line.strip('\n\r')
        posFiles.append(line + '-pos.h5')

UHDB = sys.argv[-1]

hf = h5py.File(posFiles[0], 'r')
print("open", posFiles[0], file=sys.stderr)
pos  = np.array(hf.get('pos'), dtype=int)
ch = np.array(hf.get('chr'), dtype='S5')
cnt = np.array(hf.get('cnt'), dtype=int)
refA =  np.array(hf.get('refA'), dtype='S1')
DATA = np.zeros(cnt.shape, dtype=int)
DATA += cnt
hf.close()
for fn in posFiles[1:]:
    print("fn", fn, file=sys.stderr)
    hf = h5py.File(fn, 'r')
    pos1 = np.array(hf.get('pos'), dtype=int)
    refA1 = np.array(hf.get('refA'), dtype='S1')
    assert(pos.all() == pos1.all())
    DATA += np.array(hf.get('cnt'),dtype=int)
    hf.close()

S = DATA.shape
FPOS = open(UHDB + '-pos.txt', 'w') 
FPOS.write( '\t'.join('Chromosome Position RefAllele A C G T'.split(' ')) + '\n')
for n in range(S[0]):
    FPOS.write( '\t'.join(map(str, [ch[n], pos[n], refA[n]]+list(DATA[n,:]))) + '\n')
FPOS.close()
