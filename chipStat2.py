#!/usr/bin/env python
import os, sys
from collections import defaultdict, Counter
from glob import glob
from time import time
import numpy as np

if len(sys.argv) < 2:
    print("Usage: chipStat2.py <suffix>")
    exit()

suffix = sys.argv[1]
print(suffix, file=sys.stderr)
fn = sorted(glob(suffix +"*"))
print(fn, file=sys.stderr)
F = open(fn[0], 'rb')
stats = np.load(F, allow_pickle=True)
F.close()
print('loaded', fn[0], file=sys.stderr)

for f in fn[1:]:
    with open(f, 'rb') as F:
        stats += np.load(F, allow_pickle=True)

F=open('chip-statSummary.stat', 'wb')
np.save(F, stats)
F.close()

print("Done", file=sys.stderr)
