#!/usr/bin/env python
import sys,os
import numpy as np


bimFn='small.bim'
pedFn='small.ped'
flipFn='positionsMatch-ivanCompare-flp.txt'
bimOut = 'flipped.bim'
pedOut = 'flipped.ped'

if len(sys.argv) >5:
    bimFn  = sys.argv[1]
    pedFn  = sys.argv[2]
    flipFn = sys.argv[3]
    bimOut =  sys.argv[4]
    pedOut =  sys.argv[5]
    
bim = {}
flip = {}

with open(flipFn) as f:
    for l in f:
        cs = l.strip("\n\r").split("\t")
        if cs[5] == 'Y' and cs[2] == '1.0':
            cs[5] = 'N'
        flip[(cs[0], cs[1])] = cs[5]

flipM = {i:k for i,k in zip('ACGT0', 'TGCA0')}

F = []
f1 = open(bimOut, 'w')
with open(bimFn, 'r') as f:
    for l in f:
        cs = l.strip("\n\r").split("\t")
        if (cs[0],cs[3]) in flip:
            if flip[(cs[0],cs[3])] == 'Y':
                cs[4] = flipM[cs[4]]
                cs[5] = flipM[cs[5]]
                F.append(True)
            else:
                F.append(False)
        f1.write("\t".join(cs[:6])+"\n")

f1.close()

f1 = open(pedOut, 'w')
with open(pedFn, 'r') as f:
    for l in f:
        cs = l.strip("\n\r").split("\t")
        cs[6:] = [flipM[x] if F[i/2] else x 
                  for i,x in enumerate(cs[6:])]
        f1.write("\t".join(cs)+"\n")
f1.close()
