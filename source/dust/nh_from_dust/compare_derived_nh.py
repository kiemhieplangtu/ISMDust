import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md

from restore             import restore

##================= MAIN ========================##
# Read data
cols    = ['idx', 'src', 'l', 'b', 'nhi', 'nhier', 't353', 'tauer', 'tNH', 'tNHer', 'ebv', 'ebver', 'eNH', 'eNHer', 'rad', 'rader', 'rNH', 'rNHer', 'Av', 'oh' , 'mol', 'HIobs', 'semiAtom', 'cat1', 'cat2', 'cat3', 'ok']
fmt     = ['i',    's',   'f', 'f', 'f',  'f',     'f',     'f',     'f',    'f',    'f',    'f',    'f',    'f',    'f',    'f',    'f',    'f',    'f',  'i',    'i',    'i',    'i',    'i',       'i',    'i',    'i']
nh93sc  = md.read_info_csv(cols, fmt, fname = '../../../doc/nh_from_dust.csv', skip=4, asarray=True)
sc93    = nh93sc['src']

tau     = nh93sc['t353']
tauer   = nh93sc['tauer']

ebv     = nh93sc['ebv']
ebver   = nh93sc['ebver']

tNH     = nh93sc['tNH']
tNHer   = nh93sc['tNHer']
 
eNH     = nh93sc['eNH']  
eNHer   = nh93sc['eNHer']

rNH     = nh93sc['rNH']
rNHer   = nh93sc['rNHer']

Av      = nh93sc['Av']

atom    = nh93sc['ok']
mol     = nh93sc['mol']
HIobs   = nh93sc['HIobs']
semiAtom= nh93sc['semiAtom']

# Ratios
y       = tNH/rNH
yer     = md.uncertainty_of_ratio(tNH, rNH, tNHer, rNHer)
yy      = eNH/rNH
yyer    = md.uncertainty_of_ratio(eNH, rNH, eNHer, rNHer)

# Av as X-axis
x       = Av
xer     = 3.1*ebver

# Atomic sightlines
fltr    = np.extract([atom == 1], atom)
atomSc  = np.extract([atom == 1], sc93)
atomAv  = np.extract([atom == 1], Av)
atomXer = np.extract([atom == 1], xer)

atomY   = np.extract([atom == 1], y)
atomYer = np.extract([atom == 1], yer)

atYY    = np.extract([atom == 1], yy)
atYYer  = np.extract([atom == 1], yyer)


# Molecular sightlines
molSc   = np.extract([mol == 1], sc93)
molAv   = np.extract([mol == 1], Av)
molXer  = np.extract([mol == 1], xer)

molY    = np.extract([mol == 1], y)
molYer  = np.extract([mol == 1], yer)

molYY   = np.extract([mol == 1], yy)
molYYer = np.extract([mol == 1], yyer)


# HI-observed-only sightlines
HIobsSc   = np.extract([HIobs == 1], sc93)
HIobsAv   = np.extract([HIobs == 1], Av)
HIobsXer  = np.extract([HIobs == 1], xer)

HIobsY    = np.extract([HIobs == 1], y)
HIobsYer  = np.extract([HIobs == 1], yer)

HIobsYY   = np.extract([HIobs == 1], yy)
HIobsYYer = np.extract([HIobs == 1], yyer)

# semi-Atom sightlines
smAtSc    = np.extract([semiAtom == 1], sc93)
smAtAv    = np.extract([semiAtom == 1], Av)
smAtXer   = np.extract([semiAtom == 1], xer)

smAtY     = np.extract([semiAtom == 1], y)
smAtYer   = np.extract([semiAtom == 1], yer)

smAtYY    = np.extract([semiAtom == 1], yy)
smAtYYer  = np.extract([semiAtom == 1], yyer)

###############################################
#
# Plot
#
###############################################
mpl.rcParams['axes.linewidth'] = 2
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

c1 = 'r'
c2 = 'purple'
c3 = 'b'
c4 = 'k'

f,(ax1, ax2) = plt.subplots(2, 1, gridspec_kw = {'wspace':0, 'hspace':0.02}, figsize=(16,12))
fts          = 32
labelsize    = 28
majorlght    = 12
minorlght    = 8
lgsize       = 14

major_xticks = np.arange(0., 10., 1.)
minor_xticks = np.arange(0., 10., 0.25)
major_yticks = np.arange(0., 10., 0.5)
minor_yticks = np.arange(0., 10., 0.25)

# ax1.errorbar(x, yy, xerr=xer, yerr=yyer,  color='darkgrey', marker='o', ls='None', markersize=6, markeredgecolor='darkgrey', markeredgewidth=1, capsize=0., label='') 

ax1.errorbar(atomAv, atomY, xerr=atomXer, yerr=atomYer,  color=c1, marker='o', ls='None', markersize=8, markeredgecolor=c1, markeredgewidth=1, label='') 
xerb1, = ax1.plot(atomAv, atomY, color=c1, marker='o', ls='None', markersize=10, markeredgecolor=c1, markeredgewidth=1, label='')

ax1.errorbar(smAtAv, smAtY, xerr=smAtXer, yerr=smAtYer,  color=c2, marker='o', ls='None', markersize=8, markeredgecolor=c2, markeredgewidth=1, label='') 
xerb2, = ax1.plot(smAtAv, smAtY, color=c2, marker='o', ls='None', markersize=10, markeredgecolor=c2, markeredgewidth=1, label='')

ax1.errorbar(molAv, molY, xerr=molXer, yerr=molYer,  color=c3, marker='o', ls='None', markersize=8, markeredgecolor=c3, markeredgewidth=1, label='') 
xerb3, = ax1.plot(molAv, molY, color=c3, marker='o', ls='None', markersize=10, markeredgecolor=c3, markeredgewidth=1, label='')

ax1.errorbar(HIobsAv, HIobsY, xerr=HIobsXer, yerr=HIobsYer,  color=c4, marker='o', ls='None', markersize=8, markeredgecolor=c4, markeredgewidth=1, label='') 
xerb4, = ax1.plot(HIobsAv, HIobsY, color=c4, marker='o', ls='None', markersize=10, markeredgecolor=c4, markeredgewidth=1, label='') 

ax1.set_title('', fontsize=0)
ax1.set_ylabel(r'$\mathrm{N_{H(\tau_{353})}/N_{H(\mathcal{R})}}$', fontsize=fts, fontweight='normal')

ax1.set_xticks(major_xticks)                                                       
ax1.set_xticks(minor_xticks, minor=True)
ax1.set_yticks(major_yticks)
ax1.set_yticks(minor_yticks, minor=True)

ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.set_xticks([0.1, 0.5, 1, 5])
ax1.set_yticks([0.5, 1, 2, 3, 4, 5])
ax1.get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
ax1.get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))

ax1.tick_params(axis='x', pad=0, labelbottom='off')
ax1.tick_params(axis='y', labelsize=labelsize)
ax1.tick_params(which='both', width=2)
ax1.tick_params(which='major', length=majorlght)
ax1.tick_params(which='minor', length=minorlght)

ax1.set_ylim(0.4, 3.75)
ax1.set_xlim(0.025, 6.)


# ax3  = ax1.twiny()
# ax3.set_xscale('log')
# ax3.set_yscale('log')

# ax3.set_xticks(9.4*np.array([0.1, 0.5, 1, 5]) )
# ax3.set_yticks([0.5, 1, 2, 3, 4, 5])
# ax3.get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
# ax3.get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))

# ax3.set_xlabel(r'$\mathrm{N_{H_2}\ [10^{20}cm^{-2}]}$', fontsize=fts, fontweight='normal')
# ax3.tick_params(axis='x', pad=5, labelbottom='off', labelsize=labelsize)
# ax3.tick_params(axis='y')
# ax3.tick_params(which='both', width=2)
# ax3.tick_params(which='major', length=majorlght)
# ax3.tick_params(which='minor', length=0.)

# ax3.set_xlim(0.025*9.4, 6.0*9.4)
# ax3.set_ylim(0.4, 3.75)

# ax1.tick_params(axis='x', pad=0, labelbottom='off')


#### Ax2 ###
major_xticks = np.arange(0., 10., 1.)
minor_xticks = np.arange(0., 10., 0.25)
major_yticks = np.arange(0., 10., 0.5)
minor_yticks = np.arange(0., 10., 0.25)

ax2.errorbar(x, y, xerr=xer, yerr=yer,  color='gray', marker='o', ls='None', markersize=6, markeredgecolor='gray', markeredgewidth=1, capsize=0., label='') 
xerb5, = ax2.plot(x, y, color='gray', marker='o', ls='None', markersize=6, markeredgecolor='gray', markeredgewidth=1, label='')

ax2.errorbar(atomAv, atYY, xerr=atomXer, yerr=atYYer,  color=c1, marker='o', ls='None', markersize=8, markeredgecolor=c1, markeredgewidth=1, label='')
xerb1, = ax2.plot(atomAv, atYY, color=c1, marker='o', ls='None', markersize=10, markeredgecolor=c1, markeredgewidth=1, label='')

ax2.errorbar(smAtAv, smAtYY, xerr=smAtXer, yerr=smAtYYer,  color=c2, marker='o', ls='None', markersize=8, markeredgecolor=c2, markeredgewidth=1, label='') 
xerb2, = ax2.plot(smAtAv, smAtYY, color=c2, marker='o', ls='None', markersize=10, markeredgecolor=c2, markeredgewidth=1, label='')

ax2.errorbar(molAv, molYY, xerr=molXer, yerr=molYYer,  color=c3, marker='o', ls='None', markersize=8, markeredgecolor=c3, markeredgewidth=1, label='') 
xerb3, = ax2.plot(molAv, molYY, color=c3, marker='o', ls='None', markersize=10, markeredgecolor=c3, markeredgewidth=1, label='')

ax2.errorbar(HIobsAv, HIobsYY, xerr=HIobsXer, yerr=HIobsYYer,  color=c4, marker='o', ls='None', markersize=8, markeredgecolor=c4, markeredgewidth=1, label='') 
xerb4, = ax2.plot(HIobsAv, HIobsYY, color=c4, marker='o', ls='None', markersize=10, markeredgecolor=c4, markeredgewidth=1, label='')

ax2.set_title('', fontsize=0)
ax2.set_ylabel(r'$\mathrm{N_{H(E(B-V))}/N_{H(\mathcal{R})}}$', fontsize=fts, fontweight='normal')
ax2.set_xlabel(r'$\mathrm{A_{V}\ [mag]}$', fontsize=fts, fontweight='normal')

ax2.set_xticks(major_xticks)                                                       
ax2.set_xticks(minor_xticks, minor=True)
ax2.set_yticks(major_yticks)
ax2.set_yticks(minor_yticks, minor=True)

ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.set_xticks([0.1, 0.5, 1, 5])
ax2.set_yticks([0.5, 1, 2, 3, 4, 5])
ax2.get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))
ax2.get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))

ax2.tick_params(axis='x', labelsize=labelsize, pad=5)
ax2.tick_params(axis='y', labelsize=labelsize)
ax2.tick_params(which='both', width=2)
ax2.tick_params(which='major', length=majorlght)
ax2.tick_params(which='minor', length=minorlght)

ax2.set_ylim(0.4, 3.75)
ax2.set_xlim(0.025, 6.)

axbox = ax2.get_position()
leg   = ax2.legend([xerb1, xerb2, xerb3, xerb4, xerb5], [r'$\mathrm{Purely\ atomic\ sightlines}$',\
	r'$\mathrm{Likely\ atomic\ sightlines}$', \
	r'$\mathrm{Molecular\ sightlines}$',\
	r'$\mathrm{Observed\ in\ HI\ only}$',\
	r'$\mathrm{Data\ points\ from\ upper\ panel}$' ], \
	 fontsize=lgsize, loc=(axbox.x0-0.115, axbox.y0+0.55), numpoints=1, handletextpad=-0.1)
leg.get_frame().set_linewidth(0.0)


plt.savefig('ratios_of_NH.eps', bbox_inches='tight', pad_inches=0.1, format='eps', dpi=400)
# plt.show()
## END - PLOT ##