#!/usr/bin/env python
from math import *
import numpy as np
from time import time
import sys
from glob import glob
posFiles = []
hpthFiles = []

if len(sys.argv) < 3:
    print("Usage: unionHDB.py HDB_1 .. HDB_n  UHDB")
    sys.exit(1)

if len(sys.argv)==3:
    # if there are only two argumantes the first is the directory name
    # with hdbs
    hpthFiles = sorted(glob(sys.argv[1] + '/*-hpth.txt'))
else:
    for line in sys.argv[1:-1]:
        print(line)
        line = line.strip()
        posFiles.append(line + '-pos.txt')
        hpthFiles.append(line + '-hpth.txt')

UHDB = sys.argv[-1]
st = time()

out = open(UHDB + '-hpth.txt', 'w')
out.write('\t'.join('Family childId haplothreadId Haplothread'.split(' ')) + '\n')
for f in hpthFiles:
    print(f)
    with open(f, 'r') as F:
        F.readline()
        for line in F:
            out.write(line)
        print(f, (time() - st), file=sys.stderr)
out.close()

