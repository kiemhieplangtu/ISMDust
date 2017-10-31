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
cols    = ['idx', 'src', 'l', 'b', 'nhi', 'nhier', 'tNH', 'tNHer', 'tNH2', 'tNH2er' , 'ebver', 'eNH', 'eNHer', 'eNH2', 'eNH2er', 'rNH', 'rNHer', 'rNH2', 'rNH2er', 'Av', 'oh' , 'noh', 'noher']
fmt     = ['i',    's',  'f', 'f', 'f',    'f',     'f',     'f',     'f',    'f',     'f',     'f',    'f',     'f',    'f',    'f',    'f',     'f',   'f',      'f',    'i',   'f',   'f'   ]
nh93sc  = md.read_info_csv(cols, fmt, fname = '../../../doc/nh2_vs_noh.csv', skip=4, asarray=True)

sc93    = nh93sc['src']
ebver   = nh93sc['ebver']

tNH2    = nh93sc['tNH2']
tNH2er  = nh93sc['tNH2er']
 
eNH2    = nh93sc['eNH2']  
eNH2er  = nh93sc['eNH2er']

rNH2    = nh93sc['rNH2']
rNH2er  = nh93sc['rNH2er']

Av      = nh93sc['Av']
oh      = nh93sc['oh']

noh     = nh93sc['noh']
noher   = nh93sc['noher']

# LOS with OH detection
fltr    = np.extract([oh == 1], oh)
xsc     = np.extract([oh == 1], sc93)
noh     = np.extract([oh == 1], noh)
noher   = np.extract([oh == 1], noher)

tNH2    = np.extract([oh == 1], tNH2)
tNH2er  = np.extract([oh == 1], tNH2er)

eNH2    = np.extract([oh == 1], eNH2)
eNH2er  = np.extract([oh == 1], eNH2er)

rNH2    = np.extract([oh == 1], rNH2)
rNH2er  = np.extract([oh == 1], rNH2er)

n       = len(fltr)

# X(OH) and errors in Unit of 1e-7
tXOH    = 10.*noh/tNH2
eXOH    = 10.*noh/eNH2
rXOH    = 10.*noh/rNH2

tXOHer  = 10.*md.uncertainty_of_ratio(noh, tNH2, noher, tNH2er)
eXOHer  = 10.*md.uncertainty_of_ratio(noh, eNH2, noher, eNH2er)
rXOHer  = 10.*md.uncertainty_of_ratio(noh, rNH2, noher, rNH2er)

## Nothing but not show the Negative N(H2) points
for i in range(n):
	if(tNH2[i] < 0.):
		tNH2[i] = tNH2[i]*1000.

	if(eNH2[i] < 0.):
		eNH2[i] = eNH2[i]*1000.

	if(rNH2[i] < 0.):
		rNH2[i] = rNH2[i]*1000.


################################################################
#
# Plot
#
################################################################

mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

mks = 8
fts = 24

c1  = 'r'
c2  = 'b'
c3  = '#363737'
mk1 = 'o'
mk2 = 'o'
mk3 = 'o'

# Plot X(OH) vs N(H2)
fig = plt.figure(figsize=(12,12))
ax  = fig.add_subplot(111); #ax.set_rasterized(True)


major_xticks = np.arange(0., 500., 10.)
minor_xticks = np.arange(0., 500., 10.)
major_yticks = np.arange(0.1, 50., 2.)                                              
minor_yticks = np.arange(0.1, 50., 1.)

xerb1, = plt.plot(tNH2, tXOH, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label=r'$X_{OH}\ from\ \tau_{353}$')
plt.errorbar(tNH2, tXOH, xerr=tNH2er, yerr=tXOHer, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label='data')

xerb2, = plt.plot(eNH2, eXOH, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label=r'$X_{OH}\ from\ E(B-V)$')
plt.errorbar(eNH2, eXOH, xerr=eNH2er, yerr=eXOHer, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label='data')

xerb3, = plt.plot(rNH2, rXOH, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label=r'$X_{OH}\ from\ \mathcal{R}$')
plt.errorbar(rNH2, rXOH, xerr=rNH2er, yerr=rXOHer, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label='data')

# c1 = mpl.patches.Ellipse((xt_nh2*1e20, xt_xoh), 1e20, 0.04, edgecolor='r', facecolor='none', linewidth=2)
# ax.add_artist(c1)

# c2 = mpl.patches.Ellipse((xe_nh2*1e20, xe_xoh), 0.8e20, 0.06, edgecolor='r', facecolor='none', linewidth=2)
# ax.add_artist(c2)

# c3 = mpl.patches.Ellipse((xr_nh2*1e20, xr_xoh), 0.3e20, 0.7, edgecolor='r', facecolor='none', linewidth=2)
# ax.add_artist(c3)

plt.title('', fontsize=0)
plt.xlabel(r'$\mathrm{N_{H_{2}}\ [cm^{-2}]} $', fontsize=fts, fontweight='normal')
plt.ylabel(r'$\mathrm{X_{OH} [10^{-7}]}$', fontsize=fts, fontweight='normal')
                                     
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=22, pad=7)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=12)
plt.tick_params(which='minor', length=6)
plt.grid(False)

plt.yscale('log')
plt.xscale('log')

plt.xlim(0.5, 300.)
plt.ylim(0.08, 15.)

axbox = ax.get_position()
leg   = plt.legend([xerb3, xerb2, xerb1], [r'$X_{OH}\ from\ \mathcal{R}$', r'$X_{OH}\ from\ E(B-V)$', r'$X_{OH}\ from\ \tau_{353}$'],\
    fontsize=22, loc=(axbox.x0+0.48, axbox.y0+0.73), numpoints=1, handletextpad=-0.3)
leg.get_frame().set_linewidth(0.0)

plt.tight_layout()

plt.savefig('xoh_vs_nh2.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
# plt.savefig('xoh_vs_nh2.png', bbox_inches='tight', pad_inches=0.03, format='png', dpi=100)
plt.show()




# Plot N(OH) vs N(H2)
fig = plt.figure(figsize=(12,12))
ax  = fig.add_subplot(111); #ax.set_rasterized(True)

major_xticks = np.arange(0., 500., 10.)
minor_xticks = np.arange(0., 500., 10.)
major_yticks = np.arange(0.1, 50., 2.)                                              
minor_yticks = np.arange(0.1, 50., 1.)

xerb1, = plt.plot(tNH2, noh, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label='')
plt.errorbar(tNH2, noh, xerr=tNH2er, yerr=noher, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label='')

xerb2, = plt.plot(eNH2, noh, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label='')
plt.errorbar(eNH2, noh, xerr=eNH2er, yerr=noher, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label='')

xerb3, = plt.plot(rNH2, noh, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label='')
plt.errorbar(rNH2, noh, xerr=rNH2er, yerr=noher, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label='')

plt.title('', fontsize=0)
plt.xlabel(r'$\mathrm{N_{H_{2}}\ [10^{20}cm^{-2}]} $', fontsize=fts, fontweight='normal')
plt.ylabel(r'$\mathrm{N_{OH} [10^{14}cm^{-2}]}$', fontsize=fts, fontweight='normal')
                                     
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=22, pad=7)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=12)
plt.tick_params(which='minor', length=6)
plt.grid(False)

plt.yscale('log')
plt.xscale('log')

plt.xlim(0.6, 300.)
plt.ylim(0.07, 10.)

# for i in range(n):
# 	plt.annotate('('+str(xsc[i])+')', xy=(tNH2[i], noh[i]), xycoords='data', xytext=(-50.,30.), textcoords='offset points', arrowprops=dict(arrowstyle="->"),fontsize=12)

axbox = ax.get_position()
leg   = plt.legend([xerb3, xerb2, xerb1], [r'$\mathrm{N_{H_2}\ from\ \mathcal{R}}$', r'$\mathrm{N_{H_2}\ from\ E(B-V)}$', r'$\mathrm{N_{H_2}\ from\ \tau_{353}}$'],\
    fontsize=20, loc=(axbox.x0-0.12, axbox.y0+0.75), numpoints=1, handletextpad=-0.3)
leg.get_frame().set_linewidth(0.0)

plt.tight_layout()

plt.savefig('noh_vs_nh2.eps', bbox_inches='tight', pad_inches=0.08, format='eps', dpi=600)
plt.show()






# Plot subplots
fig          = plt.figure(figsize=(14,10))
lbsize       = 22
ftsz         = 32

# Subplot 1
ax           = fig.add_subplot(121); #ax.set_rasterized(True)
major_xticks = np.arange(0.0, 0.4, 0.05)
minor_xticks = np.arange(0.0, 0.4, 0.025)
major_yticks = np.arange(0.0, 0.4, 0.05)
minor_yticks = np.arange(0.0, 0.4, 0.025)

major_xticks = np.arange(0., 500., 10.)
minor_xticks = np.arange(0., 500., 10.)
major_yticks = np.arange(0.1, 50., 2.)                                              
minor_yticks = np.arange(0.1, 50., 1.)

xerb1, = ax.plot(tNH2, noh, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label='')
ax.errorbar(tNH2, noh, xerr=tNH2er, yerr=noher, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label='')

xerb2, = ax.plot(eNH2, noh, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label='')
ax.errorbar(eNH2, noh, xerr=eNH2er, yerr=noher, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label='')

xerb3, = ax.plot(rNH2, noh, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label='')
ax.errorbar(rNH2, noh, xerr=rNH2er, yerr=noher, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label='')

ax.set_xlabel(r'$\mathrm{N_{H_{2}}\ [10^{20}cm^{-2}]} $', fontsize=fts, fontweight='normal')
ax.set_ylabel(r'$\mathrm{N_{OH} [10^{14}cm^{-2}]}$', fontsize=fts, fontweight='normal')

ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
ax.tick_params(axis='x', labelsize=18, pad=6)
ax.tick_params(axis='y', labelsize=18)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=11)
ax.tick_params(which='minor', length=5)
ax.grid(False)

ax.set_yscale('log')
ax.set_xscale('log')

ax.set_xlim(0.6, 300.)
ax.set_ylim(0.07, 10.)

axbox = ax.get_position()
leg   = plt.legend([xerb3, xerb2, xerb1], [r'$\mathrm{N_{H_2}\ from\ \mathcal{R}}$', r'$\mathrm{N_{H_2}\ from\ E(B-V)}$', r'$\mathrm{N_{H_2}\ from\ \tau_{353}}$'],\
    fontsize=16, loc=(axbox.x0-0.1, axbox.y0+0.745), numpoints=1, handletextpad=-0.3)
leg.get_frame().set_linewidth(0.0)

## Subplot 2
ax1          = fig.add_subplot(122); #ax.set_rasterized(True)
major_xticks = np.arange(0., 20., 10.)
minor_xticks = np.arange(0., 20., 10.)
major_yticks = np.arange(0.1, 20., 2.)                                              
minor_yticks = np.arange(0.1, 20., 1.)

xerb1, = ax1.plot(tNH2, tXOH, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label=r'$X_{OH}\ from\ \tau_{353}$')
ax1.errorbar(tNH2, tXOH, xerr=tNH2er, yerr=tXOHer, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label='data')

xerb2, = ax1.plot(eNH2, eXOH, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label=r'$X_{OH}\ from\ E(B-V)$')
ax1.errorbar(eNH2, eXOH, xerr=eNH2er, yerr=eXOHer, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label='data')

xerb3, = ax1.plot(rNH2, rXOH, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label=r'$X_{OH}\ from\ \mathcal{R}$')
ax1.errorbar(rNH2, rXOH, xerr=rNH2er, yerr=rXOHer, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label='data')

ax1.set_xlabel(r'$\mathrm{N_{H_{2}}\ [10^{20}cm^{-2}]} $', fontsize=fts, fontweight='normal')
ax1.set_ylabel(r'$\mathrm{X_{OH} [10^{-7}]}$', fontsize=fts, fontweight='normal')
                                     
# ax.set_yticks(major_yticks)                                                       
# ax.set_yticks(minor_yticks, minor=True)
ax1.tick_params(axis='x', labelsize=18, pad=6)
ax1.tick_params(axis='y', labelsize=18)
ax1.tick_params(which='both', width=1.5)
ax1.tick_params(which='major', length=11)
ax1.tick_params(which='minor', length=5)
ax1.grid(False)

ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.set_xlim(0.5, 300.)
ax1.set_ylim(0.08, 15.)

axbox = ax.get_position()
leg   = ax1.legend([xerb3, xerb2, xerb1], [r'$\mathrm{X_{OH}\ from\ \mathcal{R}}$', r'$\mathrm{X_{OH}\ from\ E(B-V)}$', r'$\mathrm{X_{OH}\ from\ \tau_{353}}$'],\
    fontsize=16, loc=(axbox.x0+0.39, axbox.y0+0.745), numpoints=1, handletextpad=-0.3)
leg.get_frame().set_linewidth(0.0)

fig.subplots_adjust(wspace=4)

plt.tight_layout()
# plt.savefig('OH_vs_H2.eps', bbox_inches='tight', pad_inches=0.08, format='eps', dpi=600)
plt.savefig('OH_vs_H2.png', bbox_inches='tight', pad_inches=0.08, format='png', dpi=600)
plt.show()