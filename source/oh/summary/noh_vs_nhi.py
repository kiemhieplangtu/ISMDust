import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md

from restore             import restore

##================= MAIN ========================##

## Infor for 93 src
dat   = md.read_info_93src_csv(fname = '../../../doc/observed_src.csv',asarray=True)
xsc   = dat['src']
xl    = dat['l']
xb    = dat['b']
nhi   = dat['nhi']  
nhier = dat['nhier']
xOH   = dat['oh']
tsg67 = dat['tsig67']
noh   = dat['noh']  
noher = dat['noher']

fltr   = np.extract([xOH == 1], xOH)
xsc1   = np.extract([xOH == 1], xsc)
xl1    = np.extract([xOH == 1], xl)
xb1    = np.extract([xOH == 1], xb)
nhi1   = np.extract([xOH == 1], nhi)
nhier1 = np.extract([xOH == 1], nhier)
noh1   = np.extract([xOH == 1], noh)
noher1 = np.extract([xOH == 1], noher)
n1     = len(fltr)

fltr   = np.extract([xOH == 0], xOH)
tsg67  = np.extract([xOH == 0], tsg67)
nhi2   = np.extract([xOH == 0], nhi)
nhier2 = np.extract([xOH == 0], nhier)

noh2   = 2.36634*3.0*tsg67*3.5

print n1

plt.errorbar(nhi1, noh1, yerr=noher1, xerr=nhier1, color='k', marker='h', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='')
plt.errorbar(nhi2, noh2, yerr=0.20*np.array(noh2), color='gray', marker='v', ls='None', uplims=True, \
markeredgewidth=1, markersize=9, capsize=2, markeredgecolor='gray', label='')


plt.ylabel(r'$\mathrm{N_{OH} [10^{14} cm^{-2}]}$', fontsize=20)
plt.xlabel(r'$\mathrm{N_{HI} [10^{20} cm^{-2}]}$', fontsize=20)
plt.xscale('log')
plt.yscale('log')
plt.savefig('oh_vs_nhi.png', bbox_inches='tight', pad_inches=0.03, format='png', dpi=150)
plt.show()