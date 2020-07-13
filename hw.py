#!/usr/bin/env python
#import h5py
import sys, os
import numpy as np
from numpy import *
from scipy.stats import *
from numpy.random import *

def isPseudoAutosomalX( pos, par1par2 ):
        pars = list(map(int,par1par2.split(",")))
        return (int(pos) >= pars[0]) and (int(pos) < pars[1]) or (int(pos) >= pars[2]) and (int(pos) < pars[3])


def HW(AA,Aa,aa):
    n = 1.*(AA + Aa + aa)
    if n == 0:
        return 0.0
    p = 1.*(2*AA + Aa)/(2*n)
    q = 1 - p
    EAA = n*p**2
    EAa = n*2*p*q
    Eaa = n*q**2
    return (AA-EAA)**2/EAA + (Aa - EAa)**2/EAa + (aa - Eaa)**2/Eaa

def HW1(AA,Aa,aa):
    n = 1.*(AA + Aa + aa)
    if n == 0:
        return 0.0
    p = 1.*(2*AA + Aa)/(2*n)
    q = 1 - p
    EAA = n*p**2
    EAa = n*2*p*q
    Eaa = n*q**2
    return np.array([EAA,EAa,Eaa])

def HW_X(mAA, maa, fAA, fAa, faa):
    """
    Expectation calculation for X chromosome pure
    male counts: mAA, maa
    female counts: fAA, fAa, faa
    """
    nm = mAA + maa
    nf = fAA + fAa + faa
    n = nm + nf
    if n == 0:
        return np.array([0,0,0,0,0])
    f = 1.*nm/n
    p = (2*fAA + fAa + mAA)/(2.0*nf+nm)
    q = 1-p
    emA = n*f*p
    ema = n*f*(1-p)
    efAA = n*(1-f)*p**2
    efAa = 2*n*(1-f)*p*(1-p)
    efaa = n*(1-f)*(1-p)**2
    return np.array([emA, ema, efAA, efAa, efaa])
    
def randomSampling( cnt, genF, smpl_size=10000, flagX=False ):
    s = sum(cnt)
    eCnt = s*np.array(genF)

    T = sum( [1.*(c-e)*(c-e)/e for c,e in zip(cnt,eCnt)] )

    x = multinomial( s, genF, size=smpl_size )

    w = (x - eCnt)*(x - eCnt) / (1.*eCnt)
    n = sum( sum(w,1) > T )
        
    pv = (1.*n)/smpl_size
    print ("randomSampling", ','.join( [str(x) for x in cnt] ), ','.join( ['{0:.2f}'.format( x ) for x in eCnt] ), pv, file=sys.stderr)
    return pv

def Chi2_test( cnt, eCnt ):
        df = len(cnt) - 2        

        T = sum([ (c-e)*(c-e)/e for c,e in zip(cnt,eCnt) if e != 0])
        pv = 1. - chi2.cdf( sum(T), df )
        print ("Chi2_test", ','.join( [str(x) for x in cnt] ), ','.join( ['{0:.2f}'.format( x ) for x in eCnt] ), pv, file=sys.stderr)
        return pv

def Chi2_options( cnt, eCnt, genF, X=False ):
        print("Chi2_options", file=sys.stderr)
        if sum(np.array(eCnt)<5) < 1 :
                return Chi2_test( cnt, eCnt )

        if (cnt[0] == 0 and eCnt[0] < 1) or (cnt[2] == 0 and eCnt[2] < 1 ):
                return 1.

        return randomSampling( cnt, genF, flagX=X )


def Test( cnt, eCnt, genF, X=False ):
        
    pv = Chi2_options( cnt, eCnt, genF, X )
    
    return pv

def pval_count_autosome( cnt ):
        # cnt: [RR, RA, AA]
        N = sum(cnt)
        p = (1.0*cnt[1] + 2.*cnt[2])/(2.*N)

        genF = [(1-p)*(1-p), 2.*(1-p)*p, p*p]
        eCnt = [ N*x for x in genF]

        pv = Test( cnt, eCnt, genF ) 

        return pv, eCnt

def pval_count_X( cnt ):
    # cnt: [mRR, mAA, fRR, fRA, fAA]
    nm = sum(cnt[:2]) 
    nf = sum(cnt[2:]) 
    N = nm + nf
    if N == 0:
        return 1, np.array([0,0,0,0,0])
    f = 1.*nm/N
    p = (2*cnt[2] + cnt[3] + cnt[0])/(2.0*nf+nm)

    #genF = [(1-p)*(1-p)/2.+(1-p)/2., (1-p)*p + p/2., p*p/2.]
    genF = [f*p, f*(1-p),(1-f)*p*p, 2*(1-f)*p*(1-p), (1-f)*(1-p)*(1-p)]
    eCnt = [ N*x for x in genF]

    pv = Test( cnt, eCnt, genF, True )

    return pv, eCnt

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

