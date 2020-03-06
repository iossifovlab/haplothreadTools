#!/usr/bin/env python
import sys, os

quads = []
with open('quadIds.txt', 'r') as f:
    for l in f:
        quads.append(l.strip("\n\r"))

out=open('quads.ped', 'w')

with open('joined.ped', 'r') as f:
    for l in f:
        if l[:5] in quads:
            out.write(l)

out.close()
