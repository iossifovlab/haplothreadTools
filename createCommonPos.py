#!/usr/bin/env python
from math import *
import numpy as np
import sys
from collections import defaultdict
from time import time

pos= [defaultdict(list), defaultdict(list)]

# put smaller file first in the arguments
if len(sys.argv) < 4:
    print "Usage: createCommonPos.py <file1> <file2> <common file>"
    sys.exit(1)

fns = sys.argv[1:3]
out = open(sys.argv[3], 'w')
st = time

with open(fns[0], 'r') as f:
    HD = {v:i for i,v in enumerate(f.readline().strip('\n\r').split('\t'))}
    n = 0
    for l in f:
        cs = l.strip('\n\r').split('\t')
        alleles = set([i for i,v in enumerate(cs[3:]) if int(v) > 0])
        pos[0][tuple(cs[:2])]=[n,alleles]
        n += 1

    print >>sys.stderr, "read pos file", fns[0]

with open(fns[1], 'r') as f:
    HD = {v:i for i,v in enumerate(f.readline().strip('\n\r').split('\t'))}
    n = 0
    for l in f:
        cs = l.strip('\n\r').split('\t')
        pos[1][tuple(cs[:2])]=cs
        n += 1
        
    print >>sys.stderr, "read pos file", fns[1]

# output only positions the have the same chr, pos and the same alleles
# we do not want introduce third allele
n = 0
for k,v in sorted(pos[0].items(), key=lambda x: x[1][0]):
    if k in pos[1]:
        s1 = pos[0][k][1]
        s2 = set([i for i,v in enumerate(pos[1][k][3:]) if int(v) > 0])
        if s1 < s2 or s2 < s1 or s1 == s2:
            out.write('\t'.join(k) + '\n')
    n += 1
    if n % 10000==0:
        print >>sys.stderr, "Processed %d lines" % n

out.close()
print >>sys.stderr, "Done"


