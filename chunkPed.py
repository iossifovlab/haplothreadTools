#!/usr/bin/env python
import os,sys

if len(sys.argv) <4:
    print("Usage: chunkPed.py <ped file> <fam list> <chunk number> <famId length>")
    sys.exit(1)

    
pedFn = sys.argv[1] 
famFn = sys.argv[2]
outFn = sys.argv[3]
famIdLen = int(sys.argv[4])

chunkN = int(outFn.split("/")[1].split(".")[0])

# works for sorted ped file

famList = []
with open(famFn, 'r') as f:
    for l in f:
        famList.append(l.strip('\n\r'))

famList = famList[50*chunkN:50*(chunkN+1)]

inList = False
with open(pedFn, 'r') as fin, open(outFn, 'w') as fout:
    for l in fin:
        if l[:famIdLen] in famList:
            fout.write(l)
            inList = True
        elif inList:
            break

        
