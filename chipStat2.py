#!/usr/bin/env python
import os, sys
from collections import defaultdict, Counter
from glob import glob
from time import time
import numpy as np

if len(sys.argv) < 2:
    print "Usage: chipStat2.py <suffix>"
    exit()

suffix = sys.argv[1]
print >>sys.stderr, suffix
fn = sorted(glob(suffix +"*"))
print >>sys.stderr, fn
F = open(fn[0], 'rb')
stats = np.load(F)
F.close()
print >>sys.stderr, 'loaded', fn[0]

for f in fn[1:]:
    with open(f, 'rb') as F:
        stats += np.load(F)

F=open('chip-statSummary.stat', 'wb')
np.save(F, stats)
F.close()

print >>sys.stderr, "Done"
