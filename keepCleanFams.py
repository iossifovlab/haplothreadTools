#!/usr/bin/env python
import sys, os

if len(sys.argv) < 5:
    print "Usage: keepCleanFams.py <ped in> <fam ids> <ped out> <id length>"  
    sys.exit(1)

inFn=sys.argv[1]
famFn = sys.argv[2]
outFn=sys.argv[3]
L = int(sys.argv[4])

fams = []
with open(famFn, 'r') as f:
    for l in f:
        fams.append(l.strip("\n\r"))

out=open(outFn, 'w')

with open(inFn, 'r') as f:
    for l in f:
        if l[:L] in fams:
            out.write(l)

out.close()
