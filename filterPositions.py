#!/usr/bin/env python
import os, sys
from collections import defaultdict

"""
autosome only
no indels
MAF > 0
percent called > 95%
HW >= 0.00001
for positiones that satisfy these criteria columns "percentCalled   MAF     HW_p" from statSammureT.txt and statSummaryM.txt are pasted the output.
"""
#chrom   pos     majorA  minorA  MM_N    Mm_N    mm_N    00_N    majorA_N        minorA_N        0_N      percentCalled   MAF     HW_p

fn = 'uniquePos.txt'
if len(sys.argv) <4:
    print("Usage: filterPositions.py <in unique pos file> <trio summary f> <multiplex summary f>")
    exit()
    
uFn = sys.argv[1]
tFn = sys.argv[2]
mFn = sys.argv[3]
st = True

def loadFlie(fn,st):
    res = defaultdict()
    with open(fn, 'r') as f:
        cs = f.readline().strip('\n\r').split('\t')
        if st:
            print('\t'.join(cs + cs[-3:] + cs[-3:]))
            st = False
        HD = {v:i for i,v in enumerate(cs)}
        for l in f:
            cs = l.strip('\n\r').split('\t')
            ch = int(cs[HD['chrom']])
            pos = int(cs[HD['pos']]) 
            res[(ch,pos)] = cs
    return res,st,HD

unique,st,HD = loadFlie(uFn,st)
trio,st,HD = loadFlie(tFn,st)
multiplex,st,HD = loadFlie(mFn,st)

for k in sorted(unique):
    q = unique[k]
    t = trio[k]
    m = multiplex[k]
    majorMatch = (q[HD['majorA']] == t[HD['majorA']] and q[HD['majorA']] == m[HD['majorA']])
    minorMatch = (q[HD['minorA']] == t[HD['minorA']] and q[HD['minorA']] == m[HD['minorA']])
    ch = (k[0] < 23)
    no_indels = (len(q[HD['majorA']]) == 1 and len(q[HD['minorA']]) == 1
    and len(t[HD['majorA']]) == 1 and len(t[HD['minorA']]) == 1
    and len(m[HD['majorA']]) == 1 and len(m[HD['minorA']]) == 1)
    percent_called = (float(q[HD['percentCalled']]) >= 95 and float(t[HD['percentCalled']]) >= 95
    and float(m[HD['percentCalled']]) >= 95 )
    maf = (float(q[HD['MAF']]) > 0 and float(t[HD['MAF']]) > 0 and float(m[HD['MAF']]) > 0 )
    hw = (float(q[HD['HW_p']]) > 0.00001 and float(t[HD['HW_p']]) > 0.00001 and float(m[HD['HW_p']]) > 0.00001)
    if ch and no_indels and maf and hw and percent_called and majorMatch and minorMatch:
        print('\t'.join(q + t[-3:] + m[-3:]))

            
        
