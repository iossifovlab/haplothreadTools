from pylab import *
import pylab as plt
import numpy as np
interactive(1)
from matchHT import loadHDB, compareHT
import sys

if len(sys.argv) < 2:
    print("Usage: statsHDB.py <HDB>")
    sys.exit(1) 

HDB = sys.argv[1]

HT, pos = loadHDB(HDB)

#### define list of haplothreads ####
keys = sorted(HT.keys())

mt = [k for k in keys if 'MT' in k]
ft = [k for k in keys if 'FT' in k]
mnt = [k for k in keys if 'MNT' in k]
fnt = [k for k in keys if 'FNT' in k]

### plot histograms ###

def htPlotFamily(HT,keys):
    kN = len(keys)
    nint = []
    nint0 = []
    tL = []
    tL0 = []
    ML = []
    KL = []
    K4L = []
    maxL = 0
    gl = []
    for k in range(kN):
        for l in range(kN):
            _, _, L, L0 = compareHT(HT, pos, keys[k], keys[l])
            gl.append(L)
            tL0.append("{:,}".format(sum(L0)))
            tL.append("{:,}".format(sum(L)))
            nint.append("{:,}".format(len(L)))
            nint0.append("{:,}".format(len(L0)))
            x = sum([1 for i in L if i >10**6])
            y = sum([i for i in L if i >10**6])
            ML.append([int(x), 100.*y/sum(L)])
            x = sum([1 for i in L if i >10**5])
            y = sum([i for i in L if i >10**5])
            KL.append([int(x),100.*y/sum(L)])
            x = sum([1 for i in L if i >10**4])
            y = sum([i for i in L if i >10**4])
            K4L.append([int(x),100.*y/sum(L)])
            if max(L) > maxL:
                maxL = max(L)

    bins = 50
    rang = (0,np.log10(maxL))
    font_size = 7
    plt.rcParams['xtick.labelsize'] = font_size
    plt.rcParams['ytick.labelsize'] = font_size
    n = 0
    xlimits = []
    ylimits = []
    ax = []
    fig,axis1 = plt.subplots(nrows=8,ncols=8)
    subplots_adjust(wspace=.3, hspace=.3)

    fig.set_figwidth(kN*22/8.)
    fig.set_figheight(kN*22/8.)
    for k in range(kN):
        for l in range(kN):
            plt.subplot(kN,kN, n+1)
            ax.append(plt.gca())
            plt.title(" ".join(keys[k][1:]) +  ", " + " ".join(keys[l][1:]), fontsize=font_size)
            plt.hist(np.log10(gl[n]), bins = bins, range=rang, alpha=0.6, edgecolor='none')
            text(.4,.75,'I:'+nint[n]+'\n'\
                 +'I_1Mb:'+str(ML[n][0])+'('+'%.2f' % ML[n][1]+'%)'+'\n'\
                 +'I_100Kb:'+str(KL[n][0])+'('+'%.2f' % KL[n][1]+'%)'+'\n'\
                 +'I_10Kb:'+str(K4L[n][0])+'('+'%.2f' % K4L[n][1]+'%)'+'\n'\
                 +'L:'+tL[n] \
                 ,ha='left',va='center',\
                 transform=ax[n].transAxes, fontsize=font_size)
            xlimits.append(ax[n].get_xlim())
            ylimits.append(ax[n].get_ylim())
            n += 1
    xlim = max([x[1] for x in xlimits])
    ylim = max([y[1] for y in ylimits])
    for i in range(len(ax)):
        ax[i].set_xlim([0,xlim])
        if ylimits[i][1] > 50:
            ax[i].set_ylim([0,ylim])

    savefig(HDB, bbox_inches='tight')
    plt.show()
    plt.close(fig)
    return [nint, nint0, tL, tL0, ML, KL]


nint, nint0, tL, tL0, ML, KL = htPlotFamily(HT,keys)

