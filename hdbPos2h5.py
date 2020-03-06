#!/usr/bin/env python

import sys, os, optparse, gzip
import h5py
import numpy as np
from time import time

if len(sys.argv) < 2:
    print("Usage: hdbPos2h5.py <HDB>")

HDB=sys.argv[1]

fin = HDB.replace('hpth','pos')
fout = HDB.replace('hpth.txt','pos.h5')
N = 11000000
t = time()
with open(fin, 'r') as f:
    cs = f.readline().strip('\n\r').split('\t')
    data_chr = np.zeros((N,), dtype='S5')
    data_pos = np.zeros((N,), dtype=int)
    data_refA = np.zeros((N,), dtype='S1')
    data_cnt = np.zeros((N,4), dtype=int)
    k = -1
    for l in f:
        k += 1
        cs = l.strip('\n\r').split('\t')
        data_chr[k] = cs[0]
        data_pos[k] = cs[1]
        data_refA[k] = cs[2]
        data_cnt[k] = cs[3:]
        #if k > 10:
        #    break

print(time() - t)

t = time()
k += 1
hf = h5py.File(fout, 'w')
hf.create_dataset('chr', data=data_chr[:k])
hf.create_dataset('pos', data=data_pos[:k])
hf.create_dataset('cnt', data=data_cnt[:k])
hf.create_dataset('refA', data=data_refA[:k])
hf.close()
print(time() - t)

"""
t = time()
hf = h5py.File(fout, 'r')
keys = hf.keys()
data = {k:np.array(hf.get(k),dtype=hf.get(k).dtype) for k in keys}
hf.close()
print time() - t

t = time()
DATA = np.zeros((N,4), dtype=int)
for n in range(2380):
    hf = h5py.File(fout, 'r')

    DATA += np.array(hf.get('cnt'),dtype=int)
    hf.close()
print time() - t
                
"""
