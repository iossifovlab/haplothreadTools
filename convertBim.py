#!/usr/bin/env python

import sys, os

if len(sys.argv) < 3:
    print "Usage: convertBim.py <old bim> <helper file>"
    exit()

bimFn = sys.argv[1]
print >>sys.stderr,  bimFn
helpFn = sys.argv[2]
print >>sys.stderr,  helpFn

bim = []
n = 0
with open(bimFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        bim.append(cs)
        n += 1

helper={}
M = {x:y for x,y in zip(list('ACGT0'),list('TGCA0'))}
with open(helpFn, 'r') as f:
    for l in f:
        cs = l.strip('\n\r').split('\t')
        helper[cs[2]] = cs
out = []
for i in range(len(bim)):
    snpId = bim[i][1]         
    if snpId in helper and len(helper[snpId][0]) < 6:
        out.append([helper[snpId][0],
                    snpId, bim[i][2],
                    helper[snpId][1],
                    helper[snpId][4],
                    helper[snpId][3]])
        print >>sys.stderr, bim[i], helper[snpId][4], helper[snpId][3]
        
    else:
        out.append(bim[i])
        out[i][3] = -1
        out[i][0] = 'chr'+out[i][0]
for i in range(len(out)):
    print '\t'.join(map(str, out[i]))
    



