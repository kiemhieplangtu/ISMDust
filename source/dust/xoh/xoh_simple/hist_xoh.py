import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md
import copy

from numpy               import array
from restore             import restore
from plotting            import cplot
from mpfit               import mpfit

## Read info of XOH from tau, Ebv, Radiance #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 4/2017
 # Author Van Hiep ##
def read_xoh(fname = 'xoh_from_tau.txt'):
	cols = ['idx', 'src', 'l', 'b', 'xoh', 'xoher', 'nh2', 'nh2er', 'nhi', 'nhier', 'nh', 'nher', 'cnm', 'cnmer', 'noh', 'noher', 'av', 'aver']
	fmt  = ['i',   's',    'f','f', 'f'    , 'f'  , 'f'  , 'f'     , 'f'  , 'f'    , 'f'  , 'f'  , 'f'  , 'f'   , 'f'   , 'f'   , 'f'   , 'f' ]
	data = restore(fname, 4, cols, fmt)
	dat  = data.read(asarray=True)
	xoh  = dat['xoh']

	return dat['src'], dat['l'], dat['b'], dat['xoh'], dat['xoher'], dat['nh2'], dat['nh2er'], dat['nhi'], dat['nhier'], dat['nh'], dat['nher'], \
	dat['cnm'], dat['cnmer'], dat['noh'], dat['noher'], dat['av'], dat['aver']

#================= MAIN ========================#
t_src, t_xl, t_xb, t_xoh, t_xoher, t_nh2, t_nh2er, t_nhi, t_nhier, t_nh, t_nher, t_cnm, t_cnmer, t_noh, t_noher, t_av, t_aver = read_xoh(fname = 'xoh_from_tau.txt')
e_src, e_xl, e_xb, e_xoh, e_xoher, e_nh2, e_nh2er, e_nhi, e_nhier, e_nh, e_nher, e_cnm, e_cnmer, e_noh, e_noher, e_av, e_aver = read_xoh(fname = 'xoh_from_ebv2011_plot.txt')
r_src, r_xl, r_xb, r_xoh, r_xoher, r_nh2, r_nh2er, r_nhi, r_nhier, r_nh, r_nher, r_cnm, r_cnmer, r_noh, r_noher, r_av, r_aver = read_xoh(fname = 'xoh_from_radiance_plot.txt')

print len(t_xoh)
print len(e_xoh)
print len(r_xoh)
print t_xoh


### X(OH) to e-7 ###
r_xoh   = r_xoh*10.
t_xoh   = t_xoh*10.
e_xoh   = e_xoh*10.
r_xoher = r_xoher*10.
t_xoher = t_xoher*10.
e_xoher = e_xoher*10.


filtr  = (t_src=='3C132')  ## filter
xt_av  = t_av[filtr]
xt_xoh = t_xoh[filtr]
xt_nh2 = t_nh2[filtr]

xe_xoh = e_xoh[filtr]
xe_nh2 = e_nh2[filtr]

xr_xoh = r_xoh[filtr]
xr_nh2 = r_nh2[filtr]

print '3C132'
print xt_av
print xt_xoh


## For Plotting ##
fts          = 42
labelsize    = 28
majorlght    = 9
minorlght    = 5
lgsize       = 34

min_y_data, max_y_data = np.min(r_xoh*100.), np.max(r_xoh*100.)
min_y_data, max_y_data = 0., 90.
binsize                = 6.25
num_y_bins             = np.floor((max_y_data - min_y_data) / binsize)
num_y_bins             = 13

# Axes definitions
nullfmt           = plt.NullFormatter()
left, width       = 0.1, 0.7
bottom, height    = 0.1, 0.8
bottom_h = left_h = left + width + 0.01

rect_scatter      = [left, bottom, width, height]
rect_histx        = [left, bottom_h, width, 0.4]
rect_histy        = [left_h, bottom, 0.2, height]

# Generate initial figure, scatter plot, and histogram quadrants
# start with a rectangular Figure
mpl.rcParams['axes.linewidth'] = 2.5
fig = plt.figure(1, figsize=(18, 10))

axScatter = plt.axes(rect_scatter)
axHistY = plt.axes(rect_histy)

# Remove labels from histogram edges touching scatter plot
axHistY.yaxis.set_major_formatter(nullfmt)

major_xticks = np.arange(0., 12., 1.)
minor_xticks = np.arange(0., 12., 0.25)
major_yticks = np.arange(0., 12., 1.)
minor_yticks = np.arange(0., 12., 0.5)

# Draw scatter plot
axScatter.errorbar(r_av, r_xoh, xerr=r_aver, yerr=r_xoher, color='r', marker='^', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label='$From$ $Radiance$')
xerb1, = axScatter.plot(r_av, r_xoh, color='r', marker='^', ls='None', markersize=10, markeredgecolor='r', markeredgewidth=1, label='$From$ $Radiance$')
axScatter.errorbar(t_av, t_xoh, xerr=t_aver, yerr=t_xoher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label=r'$From$ $\tau_{353}$')
xerb2, = axScatter.plot(t_av, t_xoh, color='b', marker='o', ls='None', markersize=10, markeredgecolor='b', markeredgewidth=1, label=r'$From$ $\tau_{353}$')
axScatter.errorbar(e_av, e_xoh, xerr=e_aver, yerr=e_xoher, color='k', marker='d', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='$From$ $E(B-V)$')
xerb3, = axScatter.plot(e_av, e_xoh, color='k', marker='d', ls='None', markersize=10, markeredgecolor='k', markeredgewidth=1, label='$From$ $E(B-V)$')


### data from vanDishoeck1986 ###
# yd = np.array([1.9, 1.7, 1.3, 2.2, 1.9, 1.4, 3.1, 2.5, 1.7, 5.0, 3.6, 2.4, 2.2, 0.68, 12., 3.2, 2.2, 2.2, 1.7])
# yd = yd/4.2
# xd = np.array([1.01, 0.79, 0.66, 0.95, 0.75, 0.64, 0.92, 0.73, 0.63, 0.86, 0.71, 0.62, 0.93, 2.12, 0.96, 0.8, 0.94, 0.95, 0.94])
# axScatter.plot(xd, yd, color='k', marker='x', ls='None', markersize=10, markeredgecolor='k', markeredgewidth=1, label='$From$ $E(B-V)$')
### data from vanDishoeck1986 ###

axScatter.set_xticks(major_xticks)
axScatter.set_xticks(minor_xticks, minor=True)
axScatter.set_yticks(major_yticks)
axScatter.set_yticks(minor_yticks, minor=True)
axScatter.tick_params(axis='x', labelsize=25, pad=3)
axScatter.tick_params(axis='y', labelsize=25)
axScatter.tick_params(which='both', width=2.5)
axScatter.tick_params(which='major', length=12)
axScatter.tick_params(which='minor', length=6)


axScatter.axvline(2.25, ymin=0.485, ymax=0.535, c='k', ls='-', linewidth=2)
axScatter.axvline(4.95, ymin=0.485, ymax=0.535, c='k', ls='-', linewidth=2)
axScatter.annotate(s='', xy=(4.97,4.5), xytext=(2.23,4.5), arrowprops=dict(arrowstyle='<->', linewidth=2))
axScatter.text(2.8, 4.59, r'$\mathrm{Sightlines\ with\ |b|<11^{o}}$', color='k', fontsize=32)
axScatter.text(0.8, 1., '(3C132)', color='k', fontsize=16, fontweight='bold')

print '3C132'
print xt_av
print xt_xoh

axScatter.annotate(s='', xy=(xt_av[0]-0.01, 0.+xt_xoh[0]), xytext=(1., 1.), arrowprops=dict(arrowstyle='->', linewidth=2))


axScatter.set_ylim(-0.2, 9.)
axScatter.set_xlim(0., 5.0)

axScatter.set_xlabel('$\mathrm{A_{V}}[mag]$', fontsize=36, fontweight='normal')
axScatter.set_ylabel('$\mathrm{X_{OH} [10^{-7}]}$', fontsize=36, fontweight='normal')

axbox = axScatter.get_position()
leg   = axScatter.legend([xerb1, xerb3, xerb2],\
   [r'$\mathrm{From\ \mathcal{R}}$',\
	r'$\mathrm{From}\ E(B-V)$',\
	r'$\mathrm{From\ \tau_{353}}$' ], \
	loc=(axbox.x0+0.5, axbox.y0+0.6), numpoints=1, fontsize=lgsize)
leg.get_frame().set_linewidth(0.0)

 
#### Draw y-axis histogram ####
### axHistY ###
major_xticks = np.arange(5., 20., 5.)
minor_xticks = np.arange(1., 20., 1.)
major_yticks = np.arange(0., 12., 1.)
minor_yticks = np.arange(0., 12., 0.5)

# Draw y-axis histogram
axHistY.hist(t_xoh, alpha=0.9,  label='', color='b', ls='-',  histtype='step', stacked=False, fill=False, range=(0.0,9.0), bins=13, lw=3, edgecolor='b', orientation='horizontal')
axHistY.hist(e_xoh, alpha=0.99, label='', color='k', ls='-',  histtype='step', stacked=False, fill=False, range=(0.0,9.0), bins=13, lw=3, edgecolor='k', orientation='horizontal')
axHistY.hist(r_xoh, alpha=1.0,  label='', color='r', ls='-', histtype='step',  stacked=False, fill=False, range=(0.0,9.0), bins=13, lw=3, edgecolor='r', orientation='horizontal')

axHistY.set_xlabel(r'$\mathrm{\#\ of\ sightlines}$', fontsize=36, fontweight='normal')

axHistY.set_xticks(major_xticks)                                                       
axHistY.set_xticks(minor_xticks, minor=True)
axHistY.set_yticks(major_yticks)
axHistY.set_yticks(minor_yticks, minor=True)
axHistY.tick_params(axis='x', labelsize=25, pad=3)
axHistY.tick_params(axis='y', labelsize=22)
axHistY.tick_params(which='both', width=2)
axHistY.tick_params(which='major', length=12)
axHistY.tick_params(which='minor', length=6)

axHistY.set_xlim(0., 15.)
axHistY.set_ylim(-0.2, 9.)

# plt.tight_layout()
plt.savefig('xoh_vs_av.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()

## X(OH) vs Av ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(10,10))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)

mks = 8
fts = 36

major_xticks = np.arange(0., 12., 1.)
minor_xticks = np.arange(0., 12., 0.25)
major_yticks = np.arange(0., 12., 1.)
minor_yticks = np.arange(0., 12., 0.5)

# Draw scatter plot
plt.errorbar(r_av, r_xoh, xerr=r_aver, yerr=r_xoher, color='r', marker='^', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label='$From$ $Radiance$')
xerb1, = plt.plot(r_av, r_xoh, color='r', marker='^', ls='None', markersize=10, markeredgecolor='r', markeredgewidth=1, label='$From$ $Radiance$')
plt.errorbar(t_av, t_xoh, xerr=t_aver, yerr=t_xoher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label=r'$From$ $\tau_{353}$')
xerb2, = plt.plot(t_av, t_xoh, color='b', marker='o', ls='None', markersize=10, markeredgecolor='b', markeredgewidth=1, label=r'$From$ $\tau_{353}$')
plt.errorbar(e_av, e_xoh, xerr=e_aver, yerr=e_xoher, color='k', marker='d', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='$From$ $E(B-V)$')
xerb3, = plt.plot(e_av, e_xoh, color='k', marker='d', ls='None', markersize=10, markeredgecolor='k', markeredgewidth=1, label='$From$ $E(B-V)$')

plt.title('', fontsize=0)
plt.xlabel('$\mathrm{A_{V}}[mag]$', fontsize=36, fontweight='normal')
plt.ylabel('$\mathrm{X_{OH} [10^{-7}]}$', fontsize=36, fontweight='normal')
                                     
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

plt.xlim(0.2, 5.0)
plt.ylim(0.08,20.)

# for i in range(len(t_src)):
# 	# if(oh[i] > 0):
# 	plt.annotate('('+str(t_src[i])+')', xy=(t_nh2[i], t_xoh[i]), xycoords='data',
#             xytext=(-50.,30.), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->"),fontsize=12,
#             )

axbox = ax.get_position()
leg   = plt.legend([xerb1, xerb2, xerb3], [r'$X_{OH}\ from\ \tau_{353}$', r'$X_{OH}\ from\ E(B-V)$', r'$X_{OH}\ from\ \mathcal{R}$'],\
    fontsize=14, loc=(axbox.x0+0.5, axbox.y0+0.7), numpoints=1)
leg.get_frame().set_linewidth(0.0)

plt.tight_layout()

# plt.savefig('xoh_vs_av.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.savefig('xoh_vs_av.png', bbox_inches='tight', pad_inches=0.03, format='png', dpi=100)
plt.show()
## END - PLOT ##



## X(OH) vs NH2 ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(12,12))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)

mks = 10
fts = 42

c3  = 'k'
c2  = 'b'
c1  = 'purple'
mk1 = '^'
mk2 = 'd'
mk3 = 'h'

major_xticks = np.arange(0., 500., 10.)
minor_xticks = np.arange(0., 500., 10.)
major_yticks = np.arange(0.1, 50., 2.)                                              
minor_yticks = np.arange(0.1, 50., 1.)

xerb1, = plt.plot(t_nh2*1e20, t_xoh, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label=r'$X_{OH}\ from\ \tau_{353}$')
plt.errorbar(t_nh2*1e20, t_xoh, xerr=t_nh2er*1e20, yerr=t_xoher, color=c1, marker=mk1, ls='None', markersize=mks, markeredgecolor=c1, markeredgewidth=1, label='data')

xerb2, = plt.plot(e_nh2*1e20, e_xoh, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label=r'$X_{OH}\ from\ E(B-V)$')
plt.errorbar(e_nh2*1e20, e_xoh, xerr=e_nh2er*1e20, yerr=e_xoher, color=c2, marker=mk2, ls='None', markersize=mks, markeredgecolor=c2, markeredgewidth=1, label='data')

xerb3, = plt.plot(r_nh2*1e20, r_xoh, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label=r'$X_{OH}\ from\ \mathcal{R}$')
plt.errorbar(r_nh2*1e20, r_xoh, xerr=r_nh2er*1e20, yerr=r_xoher, color=c3, marker=mk3, ls='None', markersize=mks, markeredgecolor=c3, markeredgewidth=1, label='data')

c1 = mpl.patches.Ellipse((xt_nh2*1e20, xt_xoh), 1e20, 0.04, edgecolor='r', facecolor='none', linewidth=2)
ax.add_artist(c1)

c2 = mpl.patches.Ellipse((xe_nh2*1e20, xe_xoh), 0.8e20, 0.06, edgecolor='r', facecolor='none', linewidth=2)
ax.add_artist(c2)

c3 = mpl.patches.Ellipse((xr_nh2*1e20, xr_xoh), 0.3e20, 0.7, edgecolor='r', facecolor='none', linewidth=2)
ax.add_artist(c3)

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

plt.xlim(0.25*1e20, 500.0*1e20)
plt.ylim(0.08,20.)

# for i in range(len(t_src)):
# 	# if(oh[i] > 0):
# 	plt.annotate('('+str(t_src[i])+')', xy=(t_nh2[i], t_xoh[i]), xycoords='data',
#             xytext=(-50.,30.), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->"),fontsize=12,
#             )

axbox = ax.get_position()
leg   = plt.legend([xerb3, xerb2, xerb1], [r'$X_{OH}\ from\ \mathcal{R}$', r'$X_{OH}\ from\ E(B-V)$', r'$X_{OH}\ from\ \tau_{353}$'],\
    fontsize=22, loc=(axbox.x0+0.48, axbox.y0+0.73), numpoints=1)
leg.get_frame().set_linewidth(0.0)

plt.tight_layout()

plt.savefig('xoh_vs_nh2.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
# plt.savefig('xoh_vs_nh2.png', bbox_inches='tight', pad_inches=0.03, format='png', dpi=100)
plt.show()
## END - PLOT ##







sys.exit()



## N(H2) vs Av ##
plt.errorbar(t_av, t_nh2, xerr=t_aver, yerr=t_nh2er, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(e_av, e_nh2, xerr=e_aver, yerr=e_nh2er, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(r_av, r_nh2, xerr=r_aver, yerr=r_nh2er, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='data')
plt.title('N$_{H2}$ (from Hiep) vs A$_{V}$', fontsize=30)
plt.xlabel('$A_{V}$ mag', fontsize=35)
plt.ylabel('$N_{H2}$', fontsize=35)
# plt.axhline(80., xmin=0, xmax=5)
# plt.axhline(10., xmin=0, xmax=5)
plt.grid(True)
# plt.ylim(-10., 80.)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
for i in range(len(t_src)):
	# if(oh[i] > 0):
	plt.annotate('('+str(t_src[i])+')', xy=(t_av[i], t_nh2[i]), xycoords='data',
            xytext=(-50.,30.), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"),fontsize=12,
            )
plt.show()



## X(OH) vs NHI ##
plt.errorbar(t_nhi, t_xoh, xerr=t_nhier, yerr=t_xoher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(e_nhi, e_xoh, xerr=e_nhier, yerr=e_xoher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(r_nhi, r_xoh, xerr=r_nhier, yerr=r_xoher, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='data')
plt.title('X$_{OH}$ (from Hiep) vs NHI', fontsize=30)
plt.xlabel('NHI', fontsize=35)
plt.ylabel('$X_{OH}$', fontsize=35)
plt.grid(True)
# plt.ylim(-10., 80.)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)

for i in range(len(t_src)):
	# if(oh[i] > 0):
	plt.annotate('('+str(t_src[i])+')', xy=(t_nhi[i], t_xoh[i]), xycoords='data',
            xytext=(-50.,30.), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"),fontsize=12,
            )
plt.show()


## X(OH) vs CNM ##
plt.errorbar(t_cnm, t_xoh, xerr=t_cnmer, yerr=t_xoher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(e_cnm, e_xoh, xerr=e_cnmer, yerr=e_xoher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(r_cnm, r_xoh, xerr=r_cnmer, yerr=r_xoher, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='data')
plt.title('X$_{OH}$ (from Hiep) vs CNM', fontsize=30)
plt.xlabel('CNM', fontsize=35)
plt.ylabel('$X_{OH}$', fontsize=35)
plt.grid(True)
# plt.ylim(-10., 80.)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
for i in range(len(t_src)):
	# if(oh[i] > 0):
	plt.annotate('('+str(t_src[i])+')', xy=(t_cnm[i], t_xoh[i]), xycoords='data',
            xytext=(-50.,30.), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"),fontsize=12,
            )
plt.show()

## X(OH) vs NH ##
plt.errorbar(t_nh, t_xoh, xerr=t_nher, yerr=t_xoher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(e_nh, e_xoh, xerr=e_nher, yerr=e_xoher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(r_nh, r_xoh, xerr=r_nher, yerr=r_xoher, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='data')
plt.title('X$_{OH}$ (from Hiep) vs NH', fontsize=30)
plt.xlabel('NH', fontsize=35)
plt.ylabel('$X_{OH}$', fontsize=35)
plt.grid(True)
# plt.ylim(-10., 80.)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
for i in range(len(t_src)):
	# if(oh[i] > 0):
	plt.annotate('('+str(t_src[i])+')', xy=(t_nh[i], t_xoh[i]), xycoords='data',
            xytext=(-50.,30.), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"),fontsize=12,
            )
plt.show()

## X(OH) vs Av ##
plt.errorbar(t_av, t_xoh, xerr=t_aver, yerr=t_xoher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(e_av, e_xoh, xerr=e_aver, yerr=e_xoher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.errorbar(r_av, r_xoh, xerr=r_aver, yerr=r_xoher, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='data')
plt.title('X$_{OH}$ (from Hiep) vs A$_{V}$', fontsize=30)
plt.xlabel('$A_{V}$ mag', fontsize=35)
plt.ylabel('$X_{OH}$', fontsize=35)
# plt.hline((0,5),(0.8, 10))
plt.grid(True)
plt.ylim(-10., 80.)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
for i in range(len(t_src)):
	# if(oh[i] > 0):
	plt.annotate('('+str(t_src[i])+')', xy=(t_av[i], t_xoh[i]), xycoords='data',
            xytext=(-50.,30.), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"),fontsize=12,
            )
plt.show()