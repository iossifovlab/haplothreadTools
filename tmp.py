#!/usr/bin/env python

from glob import glob

files = glob("*-stats.txt")

data = [[],[],[],[],[]]

expectedCodes = "ACGTcdouMRWSYK"
expectedCodesSet = set(expectedCodes)

for i,fn in enumerate(files):
    with open(fn, 'r') as f:
        f.readline()
        for l in f:
            cs = map(int,l.strip("\n\r").split("\t")[3:])
            data[i].append(cs)
        data[i] = array(data[i])


