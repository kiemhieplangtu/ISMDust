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

# Read data
cols    = ['idx', 'src', 'l', 'b', 'nhi', 'nhier', 'tau', 'tauer', 'tNH', 'tNHer', 'ebv', 'ebver', 'eNH', 'eNHer', 'rad', 'rader', 'rNH', 'rNHer', 'Av', 'oh' , 'mol', 'obsHI', 'semiA', 'noMol1', 'noMol2', 'lowNHI', 'atomic']
fmt     = ['i',    's',  'f', 'f', 'f',    'f',     'f',     'f',    'f',   'f',    'f',    'f',    'f',     'f',    'f',    'f',    'f',    'f',   'f',  'i',   'f',   'f',      'f',     'f',      'f',      'f',      'f'   ]
nh93sc  = md.read_info_csv(cols, fmt, fname = '../../../doc/nh_from_dust.csv', skip=4, asarray=True)

sc93    = nh93sc['src']
oh      = nh93sc['oh']

print ''
rad     = nh93sc['rad']
rader   = nh93sc['rader']

rNH     = nh93sc['rNH']
rNHer   = nh93sc['rNHer']

print ''
tau     = nh93sc['tau']
tauer   = nh93sc['tauer']

tNH     = nh93sc['tNH']
tNHer   = nh93sc['tNHer']

print ''
ebv     = nh93sc['ebv']
ebver   = nh93sc['ebver']

eNH     = nh93sc['eNH']
eNHer   = nh93sc['eNHer']

print len(oh)

a, aer  = [1.39249e6, 4.9102e4] #6.6e-27, 0.66e-26, lowNHI
plt.errorbar(tau, tNH, xerr=tauer, yerr=tNHer, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='')
plt.plot([0., 500.], [0.*a, 500.*a*1e-6], 'k-')
plt.show()

a, aer  = [113.912, 4.18695]
plt.plot([0., 2.], [0.*a, 2.*a], 'k-')
plt.errorbar(ebv, eNH, xerr=ebver, yerr=eNHer, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='')
plt.show()

a, aer = [4.65435e11, 0.170979e11]   ## N(H) = a.e31.R + b, NH/e20 = a.e11.R + b
plt.plot([0., 50.], [0.*a, 50.*a*1e-11], 'k-')
plt.errorbar(rad, rNH, xerr=rader, yerr=rNHer, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='')
plt.show()