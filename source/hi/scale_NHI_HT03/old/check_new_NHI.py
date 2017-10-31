import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md

from scipy.io.idl        import readsav
from restore             import restore

### 79 sources ###
ms79sc   = md.read_info_ms_79sc(fname = '../result/nhi_lb_79src_HT03.txt', asarray=True)
sc79     = ms79sc['src']
xl79     = ms79sc['l']
xb79     = ms79sc['b']
nhi03    = ms79sc['nhi']
nhier03  = ms79sc['nhi_er']
thin03   = ms79sc['thin']
thin03er = ms79sc['thin_er']
cnm03    = ms79sc['cnm']
cnm03er  = ms79sc['cnm_er']
wnm03    = ms79sc['wnm']
wnm03er  = ms79sc['wnm_er']

new      = md.read_info_ms_79sc(fname = 'nhi_79src_recal.txt', asarray=True)

nhi      = new['nhi']
nhier    = new['nhi_er']

xfit, yfit, mu, sig, m, ea, b, eb = md.do_linODRfit(nhi03, nhi, nhier03, nhier, lguess=[1.25, 0.5], plot=False)

### Plot ###
plt.figure(1, figsize=(16, 14))

plt.plot(xfit, mu, '-r', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

plt.errorbar(nhi03, nhi, xerr=nhier03, yerr=nhier, color='r', marker='o', ls='', markersize=6, markeredgecolor='r', markeredgewidth=1, label='data')

plt.plot(xfit, mu, '-k', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu-sig, mu+sig, color='0.5', alpha=0.5)

plt.plot([0,160],[0,160], 'k--', label='$x=y$')
plt.xlabel('N(HI) - HT03')
plt.ylabel('New - N(HI)')
plt.show()