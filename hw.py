#!/usr/bin/env python
#import h5py
import sys, os
import numpy as np

def HW(AA,Aa,aa):
    n = 1.*(AA + Aa + aa)
    if n == 0:
        return 0.0
    p = (2*AA + Aa)/(2*n)
    q = 1 - p
    EAA = n*p**2
    EAa = n*2*p*q
    Eaa = n*q**2
    return (AA-EAA)**2/EAA + (Aa - EAa)**2/EAa + (aa - Eaa)**2/Eaa

def HW1(AA,Aa,aa):
    n = 1.*(AA + Aa + aa)
    if n == 0:
        return 0.0
    p = (2*AA + Aa)/(2*n)
    q = 1 - p
    EAA = n*p**2
    EAa = n*2*p*q
    Eaa = n*q**2
    return np.array([EAA,EAa,Eaa])
    
def HWM(data):
    AA = np.sum(data == 0, axis=0).astype(float)
    Aa = np.sum(data == 1, axis=0).astype(float)
    aa = np.sum(data == 2, axis=0).astype(float)
    n = (AA + Aa + aa)
    p = np.divide(AA*2 + Aa, n*2)
    q = 1 - p
    EAA = np.multiply(n, np.multiply(p,p))
    EAa = np.multiply(n*2, np.multiply(p,q))
    Eaa = np.multiply(n, np.multiply(q,q))
    res = np.divide(np.multiply(AA-EAA,AA-EAA), EAA)
    res += np.divide(np.multiply(Aa-EAa,Aa-EAa), EAa)
    res += np.divide(np.multiply(aa-Eaa,aa-Eaa), Eaa)
    return res

def getCnt(data):
    AA = np.sum(data == 0, axis=0)
    Aa = np.sum(data == 1, axis=0)
    aa = np.sum(data == 2, axis=0)
    return [AA, Aa, aa]

"""
f = h5py.File('HT38Prb/AgreAllParents-chr22.hdf5', 'r')
print >>sys.stderr, "Keys: %s" % f.keys()
a_group_key = list(f.keys())[0]
agrep = np.array(f[a_group_key],dtype='|S2').astype(float)
"""

