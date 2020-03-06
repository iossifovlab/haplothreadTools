#!/usr/bin/env python
from math import *
import numpy as np
import sys
from collections import defaultdict
from time import time

# ht chrom parts are files produced by "cut -f 6- hdb|tail -n +2"

if len(sys.argv) < 3:
    print("Usage: pasteChr.py <ht parts of chroms> <hdb out>")
    sys.exit(1)

print([sys.argv[i] for i in range(1,len(sys.argv)-1)])

files = [open(sys.argv[i], 'r') for i  in range(1,len(sys.argv)-1)]
print(sys.argv[-1])

out = open(sys.argv[-1], 'w')

data={i:defaultdict() for i in range(0,22)}

for i,f in enumerate(files):
    hd = f.readline()
    for l in f:
        cs = l.strip("\n\r").split("\t")
        data[i][tuple(cs[:3])] = cs[3]

out.write(hd)

for k in sorted(data[0].keys()):
    ht = data[0][k]
    for i in range(1,22):
        if k in data[i]:
            ht += data[i][k]
    out.write("\t".join(list(k) + [ht]) + '\n')

out.close()
for f in files:
    f.close()
    
print("Done", file=sys.stderr)
