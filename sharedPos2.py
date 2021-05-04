#!/usr/bin/env python

import os,sys,optparse
import numpy as np
from collections import defaultdict
from glob import glob

#
# Chip positions contain all chromosome, but SSC and AGRE WGS haplothreads
# done by chromosomes, therefore we first find subsets pf positions common to the chip
# and each chromosome separately.
# Following conditions must be sarisfied:
# a) In all three sets position should be biallelic;
# b) refA in SSC and AGRE should be either major_A or minor_A in the chip;
# c) position on the chip should be SNP (follows from b);
# d) mosition on the chip should be genotyped at least for 95% individuals
# e) MAF on chip should be >0
# f) HW on chip should be >0.00001
#

if len(sys.argv) <3:
    print("Usage: sharedPos2.txt <chrom number> <chip filteredPos.txt")
    exit()
    
ch = sys.argv[1]
k =int( ch[3:])
statsFn = sys.argv[2]
sscDir = '/mnt/wigclust8/data/safe/autism/HT/SSC2380/'
agreDir = '/mnt/wigclust8/data/safe/autism/HT/AGRE_WG38_AGRE859/'
sscPosFns = {k:glob(sscDir + 'SSC2380-chr' + str(k) + '-pos.txt') for k in range(1,23)}
agrePosFns = {k:glob(agreDir + 'AGRE859-chr' + str(k) + '-pos.txt') for k in range(1,23)}

#chrom   pos     majorA  minorA  MM_N    Mm_N    mm_N    00_N    majorA_N        minorA_N        0_N     percentCalled   MAFHW_p


sharedMap = {}
with open(statsFn, 'r') as f:
    HDSpark = {v:i for i,v in enumerate(f.readline().strip('\n\r').split('\t'))}
    for l in f:
        cs = l.strip('\n\r').split('\t')
        #spark[(cs[0], cs[1])] = {k:cs[i] for i,k in enumerate(HDSpark)}
        sharedMap[(cs[0], cs[1])]=cs

#Chromosome      Position        RefAllele       A       C       G       T
ALLS = 'ACGT'

ssc = {}
with open(sscPosFns[k][0], 'r') as f:
    f.readline()
    for l in f:
        cs = l.strip('\n\r').split('\t')
        # because chrom has 'chr' at the beginning
        ch = cs[0][3:]
        pos = cs[1]
        refA = cs[2]
        ssc[pos] = [cs[2]] + list(map(int,cs[3:]))

agre = {}
shared = []
with open(agrePosFns[k][0], 'r') as f:
    f.readline()
    for l in f:
        cs = l.strip('\n\r').split('\t')
        # because chrom has 'chr' at the beginning            
        ch = cs[0][3:]
        pos = cs[1]
        refA = cs[2]
        key = (ch, pos)
        agre[pos] = [cs[2]] + list(map(int, cs[3:]))
        # refA should be either major or minor allele in spark
        if (key in sharedMap and
            pos in ssc and
            (refA == sharedMap[key][HDSpark['majorA']] or
             refA == sharedMap[key][HDSpark['minorA']]) and
            refA == ssc[pos][0]):
            
            sscCNT = sorted([[v,i] for i,v in enumerate(ssc[pos][1:])],
                            key=lambda x: x[0], reverse=True)
            agreCNT = sorted([[v,i] for i,v in enumerate(agre[pos][1:])],
                             key=lambda x: x[0], reverse=True)
            # condition eliminates SSC and AGRE positions with more than two alleles
            if sscCNT[2][0] == 0 and agreCNT[2][0] == 0:
                shared.append(pos)

for p in sorted(map(int,shared)):
    #print('\t'.join(map(str,[k, p])))
    print('\t'.join(sharedMap[tuple(map(str,[k,p]))][:4]))


