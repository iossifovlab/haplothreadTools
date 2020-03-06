#!/usr/bin/env python
import sys, os
from collections import defaultdict
if len(sys.argv) < 3:
    print "Usage: convertFam.py <old fam> <helper file>"
    exit()

famFn = sys.argv[1]
print >>sys.stderr,  famFn
helpFn = sys.argv[2]
print >>sys.stderr,  helpFn

fam = []
n = 0
with open(famFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split(' ')
        fam.append(cs)
        n += 1

#family  individual      array
#11000   11000.fa        4262850452_A
#11000   11000.mo        4262850453_A
helper = {}

helper=defaultdict()
with open(helpFn, 'r') as f:
    f.readline()
    for l in f:
        cs = l.strip('\n\r').split('\t')
        helper[cs[2]]=cs[1]
        n += 1

for i in range(len(fam)):
    out=fam[i]
    for k in range(1,4):
       out[k] = helper[out[k]]
    print '\t'.join(out)
    
