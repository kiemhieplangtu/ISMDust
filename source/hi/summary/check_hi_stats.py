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

## Origin infor for 93 src
cols    = ['idx', 'src', 'l', 'b', 'nhi', 'nhier', 'thin', 'thiner', 'cnm', 'cnmer', 'wnm', 'wnmer']
fmt     = ['i',    's',   'f', 'f', 'f',  'f',     'f',     'f',     'f',    'f',    'f',    'f']
hi93sc  = md.read_info_csv(cols, fmt, fname = '../../../doc/origin_hi_93src.csv', skip=4, asarray=False)
sc93    = hi93sc['src']
nhi93   = hi93sc['nhi']
nhier93 = hi93sc['nhier']
xl93    = hi93sc['l']
xb93    = hi93sc['b']
thin93  = hi93sc['thin']
thinr93 = hi93sc['thiner']
cnm93   = hi93sc['cnm']
cnmer93 = hi93sc['cnmer']
wnm93   = hi93sc['wnm']
wnmer93 = hi93sc['wnmer']

## Old infor for 93 src
cols    = ['idx', 'src', 'l', 'b', 'nhi', 'nhier', 'thin', 'thiner', 'cnm', 'cnmer', 'wnm', 'wnmer']
fmt     = ['i',    's',   'f', 'f', 'f',  'f',     'f',     'f',     'f',    'f',    'f',    'f']
dat93   = md.read_info_csv(cols, fmt, fname = '../../../doc/old_hi_93src.csv', skip=4, asarray=False)
scx     = dat93['src']
nhix    = dat93['nhi']
nhierx  = dat93['nhier']
xlx     = dat93['l']
xbx     = dat93['b']
thinx   = dat93['thin']
thinrx  = dat93['thiner']
cnmx    = dat93['cnm']
cnmerx  = dat93['cnmer']
wnmx    = dat93['wnm']
wnmerx  = dat93['wnmer']

xdat    = np.zeros(93)
xerr    = np.zeros(93)
ydat    = np.zeros(93)
yerr    = np.zeros(93)
for i in range(93):
	idx         = sc93.index(scx[i])
	xdat[i]     = nhi93[idx]
	xerr[i]     = nhier93[idx]
	ydat[i]     = nhix[i]
	yerr[i]     = nhierx[i]
	print scx[i]==sc93[idx], nhix[i], nhi93[idx]

# plot 
plt.plot([0.,140.],[0., 140.], 'k-')        # for SPONGE
plt.plot([0.,140.],[0., 1.26*140.], 'b-')   # For MS
plt.errorbar(xdat, ydat, xerr=xerr, yerr=yerr, color='r', marker='o', ls='', markersize=6, markeredgecolor='r', markeredgewidth=1, label='data')
plt.show()

sys.exit()


## Infor for 93 src
dat   = md.read_info_93src_csv(fname = '../../../doc/observed_src.csv',asarray=True)
xsc   = dat['src']
xl    = dat['l']
xb    = dat['b']
nhi   = dat['nhi']   ## Already in ascending order
nhier = dat['nhier']
xOK   = dat['ok']
tsg67 = dat['tsig67']
n     = len(xsc)

xnhi93   = np.zeros(n)
xnhier93 = np.zeros(n)
for i in range(n):
	idx         = sc93.index(xsc[i])
	# idx         = scx.index(xsc[i])	
	xnhi93[i]   = nhi93[idx]
	xnhier93[i] = nhier93[idx]

	print sc93[idx], xl93[idx], xb93[idx], nhi93[idx], nhier93[idx], thin93[idx], thinr93[idx], cnm93[idx], cnmer93[idx], wnm93[idx], wnmer93[idx]
	# print scx[idx], xlx[idx], xbx[idx], nhix[idx], nhierx[idx], thinx[idx], thinrx[idx], cnmx[idx], cnmerx[idx], wnmx[idx], wnmerx[idx]


print xnhi93.shape
print nhi.shape

# plot 
plt.plot([0.,140.],[0., 140.], 'k-')        # for SPONGE
plt.plot([0.,140.],[0., 1.28*140.], 'b-')   # For MS

plt.errorbar(xnhi93, nhi, xerr=xnhier93, yerr=nhier, color='r', marker='o', ls='', markersize=6, markeredgecolor='r', markeredgewidth=1, label='data')
plt.show()