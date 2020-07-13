#!/usr/bin/env python
import sys, os

SPARK = "SPARK" in os.environ["PROJECT_DIR"]
if SPARK:
    famN = 9
else:
    famN = 5
quads = []
with open('quadIds.txt', 'r') as f:
    for l in f:
        quads.append(l.strip("\n\r"))

out=open('quads.ped', 'w')

with open('joined.ped', 'r') as f:
    for l in f:
        if l[:famN] in quads:
            out.write(l)

out.close()
