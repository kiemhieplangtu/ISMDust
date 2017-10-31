import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot   as plt
import matplotlib          as mpl
import numpy               as np
import module              as md
import matplotlib.gridspec as gridspec

from numpy               import array
from restore             import restore



## Read NH obtained from 3 proxies for 35 los #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 08/2017
 # Author Van Hiep ##
def read_nh_from_3proxies_94src(fname = 'nh_from_tau353.txt'):
	cols = ['src', 'l', 'b', 'prx', 'prx_er', 'nh', 'nher','nhi', 'nhier', 'nh2', 'nh2er']
	fmt  = ['s',   'f', 'f',  'f', 'f'      , 'f'   , 'f'  , 'f'  , 'f'   , 'f'  , 'f'   ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=True)
	return dat

## Read infor of 94 src, atomic or molecular, CO? OH? #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 09/2017
 # Author Van Hiep ##
def read_94src_atomic_info(fname = '94src_atomic_or_molecular.txt'):
	cols = ['src', 'l', 'b', 'atom', 'coyn', 'ohyn' ]
	fmt  = ['s',   'f', 'f',  'i',     'i',   'i'   ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=True)
	return dat['src'], dat['atom'], dat['coyn'], dat['ohyn']

###================= MAIN ========================####
## Read Infor for 93 src
dat   = md.read_info_93src_csv(fname = '../../../doc/observed_src.csv',asarray=True)
xsc   = dat['src']
xl    = dat['l']
xb    = dat['b']
nhi   = dat['nhi']   ## Already in ascending order
nhier = dat['nhier']
xOK   = dat['ok']
tsg67 = dat['tsig67']

fltr  = np.extract([xOK == 1], xOK)
xsc   = np.extract([xOK == 1], xsc)
xl    = np.extract([xOK == 1], xl)
xb    = np.extract([xOK == 1], xb)
nhi   = np.extract([xOK == 1], nhi)
nhier = np.extract([xOK == 1], nhier)
tsg67 = np.extract([xOK == 1], tsg67)
n     = len(fltr)














xsc, \
atom, \
coyn, \
ohyn    = read_94src_atomic_info(fname = '94src_atomic_or_molecular.txt')

tauInfo = read_nh_from_3proxies_94src(fname = 'nh_from_tau_94src.txt')
ebvInfo = read_nh_from_3proxies_94src(fname = 'nh_from_ebv2011_94src.txt')
radInfo = read_nh_from_3proxies_94src(fname = 'nh_from_rad_94src.txt')

info94  = md.read_nhi_93src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_scaled.txt')
cnm     = info94['cnm'] 
cnmer   = info94['cnm_er']
src94   = info94['src']  

Ebv, Ebver, Av, Src = md.read_ebv_av(fname = '../ebv2nh/data/ebv_sfd98_sf2011_for_94src.txt', sfd98=False)

nhTau   = tauInfo['nh']
nherTau = tauInfo['nher']
nhEbv   = ebvInfo['nh']
nherEbv = ebvInfo['nher']
nhRad   = radInfo['nh']
nherRad = radInfo['nher']

scTau   = tauInfo['src']
scEbv   = ebvInfo['src']
scRad   = radInfo['src']

# for i in range(len(scEbv)):
# 	string = '{:10s} {:08.4f}   {:08.4f}   {}'\
# 	.format(scEbv[i], ebvInfo['l'][i], ebvInfo['b'][i], 0)
# 	print string

for i in range(94):
	print scTau[i] == scEbv[i] == scRad[i] == Src[i] == src94[i]

indx    = np.arange(94)
nhi     = radInfo['nhi']
nhier   = radInfo['nhier']

y       = nhTau/nhRad
yer     = md.uncertainty_of_ratio(nhTau, nhRad, nherTau, nherRad)
yy      = nhEbv/nhRad
yyer    = md.uncertainty_of_ratio(nhEbv, nhRad, nherEbv, nherRad)

# x       = nhi
# xer     = nhier

x       = Av
xer     = 3.1*Ebver

# x       = indx
# xer     = indx*0.

# x       = cnm
# xer     = cnmer


## For Plotting ##

filtr     = (atom==1)  ## filter
x_atom    = x[filtr]
y_atom    = y[filtr]

xer_atom  = xer[filtr]
yer_atom  = yer[filtr]

yy_atom   = yy[filtr]
yyer_atom = yyer[filtr]

filtr     = (atom==0) & ((coyn==1) | (ohyn==1))  ## filter
x_mol     = x[filtr]
y_mol     = y[filtr]

xer_mol   = xer[filtr]
yer_mol   = yer[filtr]

yy_mol    = yy[filtr]
yyer_mol  = yyer[filtr]

print len(x_mol)
print len(yy_mol)



mpl.rcParams['axes.linewidth'] = 4
f,(ax1, ax2) = plt.subplots(2, 1, gridspec_kw = {'wspace':0, 'hspace':0.02}, figsize=(24,18))
fts          = 52
labelsize    = 32
majorlght    = 16
minorlght    = 10
lgsize       = 40

major_xticks = np.arange(0., 10., 1.)
minor_xticks = np.arange(0., 10., 0.25)
major_yticks = np.arange(0., 10., 0.5)
minor_yticks = np.arange(0., 10., 0.25)

# ax1.errorbar(x, yy, xerr=xer, yerr=yyer,  color='darkgrey', marker='o', ls='None', markersize=6, markeredgecolor='darkgrey', markeredgewidth=1, capsize=0., label='') 

ax1.errorbar(x, y, xerr=xer, yerr=yer,  color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='') 
xerb1, = ax1.plot(x, y, color='k', marker='o', ls='None', markersize=10, markeredgecolor='k', markeredgewidth=1, label='') 

ax1.errorbar(x_atom, y_atom, xerr=xer_atom, yerr=yer_atom,  color='r', marker='o', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label='') 
xerb2, = ax1.plot(x_atom, y_atom, color='r', marker='o', ls='None', markersize=10, markeredgecolor='r', markeredgewidth=1, label='')

ax1.errorbar(x_mol, y_mol, xerr=xer_mol, yerr=yer_mol,  color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='') 
xerb3, = ax1.plot(x_mol, y_mol, color='b', marker='o', ls='None', markersize=10, markeredgecolor='b', markeredgewidth=1, label='') 



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
ax1.tick_params(which='both', width=3)
ax1.tick_params(which='major', length=majorlght)
ax1.tick_params(which='minor', length=minorlght)

ax1.set_ylim(0.5, 4.5)
ax1.set_xlim(0.025, 6.)

axbox = ax1.get_position()
leg   = ax1.legend([xerb2, xerb3, xerb1], [r'$\mathrm{Atomic\ sightlines}$',\
	r'$\mathrm{Molecular\ sightlines}$',\
	r'$\mathrm{Observed\ in\ HI\ only}$' ], \
	 fontsize=lgsize, loc=(axbox.x0-0.12, axbox.y0+0.09), numpoints=1)
leg.get_frame().set_linewidth(0.0)



#### Ax2 ###
major_xticks = np.arange(0., 10., 1.)
minor_xticks = np.arange(0., 10., 0.25)
major_yticks = np.arange(0., 10., 0.5)
minor_yticks = np.arange(0., 10., 0.25)

ax2.errorbar(x, y, xerr=xer, yerr=yer,  color='gray', marker='o', ls='None', markersize=6, markeredgecolor='gray', markeredgewidth=1, capsize=0., label='') 

ax2.errorbar(x, yy, xerr=xer, yerr=yyer,  color='k', marker='h', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='$\mathrm{X=E(B-V)}$') 
xerb1, = ax2.plot(x, yy, color='k', marker='h', ls='None', markersize=10, markeredgecolor='k', markeredgewidth=1, label='$\mathrm{X=E(B-V)}$') 
ax2.errorbar(x_atom, yy_atom, xerr=xer_atom, yerr=yyer_atom,  color='r', marker='h', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label=r'$\mathrm{X=\tau_{353}}$') 
xerb2, = ax2.plot(x_atom, yy_atom, color='r', marker='h', ls='None', markersize=10, markeredgecolor='r', markeredgewidth=1, label=r'$\mathrm{X=\tau_{353}}$') 

ax2.errorbar(x_mol, yy_mol, xerr=xer_mol, yerr=yyer_mol,  color='b', marker='h', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='') 
xerb3, = ax2.plot(x_mol, yy_mol, color='b', marker='h', ls='None', markersize=10, markeredgecolor='b', markeredgewidth=1, label='') 

ax2.set_title('', fontsize=0)
ax2.set_ylabel(r'$\mathrm{N_{H(E(B-V))}/N_{H(\mathcal{R})}}$', fontsize=fts, fontweight='normal')
ax2.set_xlabel(r'$\mathrm{A_{V}[mag]}$', fontsize=fts, fontweight='normal')

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
ax2.tick_params(which='both', width=3)
ax2.tick_params(which='major', length=majorlght)
ax2.tick_params(which='minor', length=minorlght)

ax2.set_ylim(0.5, 4.5)
ax2.set_xlim(0.025, 6.)

# axbox = ax2.get_position()
# leg   = ax2.legend([xerb2, xerb3, xerb1], [r'$\mathrm{Atomic\ sightlines}$',\
# 	r'$\mathrm{Molecular\ sightlines}$',\
# 	r'$\mathrm{Observed\ in\ HI\ only}$' ], \
# 	 fontsize=23, loc=(axbox.x0-0.12, axbox.y0+0.6), numpoints=1)
# leg.get_frame().set_linewidth(0.0)


plt.savefig('ratios_of_NH.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
# plt.show()
## END - PLOT ##


sys.exit()







# Draw scatter plot
# , 'r*', mew=0, linewidth=0, linestyle='', marker='*', markerfacecolor='r', markersize=8, label='N(H) from Tau353'
axScatter.errorbar(x, y, xerr=xer, yerr=yer,  color='dimgray', marker='H', ls='None', markersize=5, markeredgecolor='dimgray', markeredgewidth=1, label=r'$\mathrm{N_{H}}$ from $\mathrm{\tau_{353}}$') 
axScatter.errorbar(x, yy, xerr=xer, yerr=yyer,  color='k', marker='^', ls='None', markersize=5, markeredgecolor='k', markeredgewidth=1, label='$\mathrm{N_{H}}$ from E(B-V)') 
# axScatter.scatter(x, y/y,  marker='d', color = 'b', edgecolor='none', s=10, alpha=1)
# axScatter.scatter(x, yy, marker='s', color = 'k', edgecolor='none', s=10, alpha=1)
# axScatter.legend(loc='upper left', fontsize=18)



axScatter.axhline(y=1., xmin=-1., xmax=10., c='k', ls='--', linewidth=2)

axScatter.set_ylim(0.0, 4.5)
axScatter.set_xlim(0., 5.5)
axScatter.grid(False)
axScatter.legend(loc='upper left', fontsize=18)

# axHistY   = fig.add_subplot(122, position=rect_histy)
major_xticks = np.arange(0., 30., 5.)
minor_xticks = np.arange(0., 30., 1.)
major_yticks = np.arange(0., 10., 1.)
minor_yticks = np.arange(0., 10., 0.25)

axScatter.set_xlabel('$\mathbf{A_{V}[mag]}$', fontsize=24, fontweight='bold')
axScatter.set_ylabel('$\mathbf{N_{H(X)}/N_{H(R)}}$', fontsize=24, fontweight='bold')

# Draw y-axis histogram
axHistY.hist(y, num_y_bins, ec='darkgrey', histtype='stepfilled', stacked=False, fill=False, lw=2, orientation='horizontal')
axHistY.hist(yy, num_y_bins, ec='k', histtype='stepfilled', stacked=False, fill=False, lw=2, orientation='horizontal')

# hist(10.*r, alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(xmin,xmax), bins=10, lw=2)

plt.grid(False)
axHistY.set_xlabel('Histogram', fontsize=18, fontweight='bold')

axHistY.set_xticks(major_xticks)                                                       
axHistY.set_xticks(minor_xticks, minor=True)
axHistY.set_yticks(major_yticks)
axHistY.set_yticks(minor_yticks, minor=True)
axHistY.tick_params(axis='x', labelsize=20, pad=3)
axHistY.tick_params(axis='y', labelsize=20)
axHistY.tick_params(which='both', width=2)
axHistY.tick_params(which='major', length=6)
axHistY.tick_params(which='minor', length=3)

axHistY.set_xlim(0, 30.)
axHistY.set_ylim(0.0, 4.5)

# plt.tight_layout()
plt.savefig('ratios_of_NH.png', bbox_inches='tight' , format='png')
plt.savefig('ratios_of_NH.eps', bbox_inches='tight' , format='eps', dpi=600)
plt.show()


sys.exit()

















## Ratio vs NHI ##
x    = nhi
xer  = nhier

## For Plotting ##
min_y_data, max_y_data = np.min(y), np.max(y)
binsize                = 0.25
num_y_bins             = np.floor((max_y_data - min_y_data) / binsize)

# Axes definitions
nullfmt           = plt.NullFormatter()
left, width       = 0.1, 0.5
bottom, height    = 0.1, 0.8
bottom_h = left_h = left + width + 0.02

rect_scatter      = [left, bottom, width, height]
rect_histx        = [left, bottom_h, width, 0.4]
rect_histy        = [left_h, bottom, 0.3, height]

# Generate initial figure, scatter plot, and histogram quadrants
fig       = plt.figure(221, figsize=(16., 10.))

axScatter = fig.add_subplot(121, position=rect_scatter)
axScatter.set_xlabel('$N_{HI} [10^{20} cm^{-2}]$', fontsize=35, fontweight='bold')
axScatter.set_ylabel('$R=N_{H-fromX}/N_{H-fromRad}$', fontsize=35, fontweight='bold')
# axScatter.set_xlim(0., 130.)
# axScatter.set_ylim(0.5, 4.)
axScatter.tick_params(axis='x', labelsize=20)
axScatter.tick_params(axis='y', labelsize=15)
axScatter.set_ylim(0.5, 4.)
axScatter.grid(True)

axHistY   = fig.add_subplot(122, position=rect_histy)
axHistY.set_xlim(0, 30.)
axHistY.set_ylim(0.5, 4.)
axHistY.set_xlabel('$Histogram$', fontsize=35, fontweight='bold')

# Remove labels from histogram edges touching scatter plot
# axHistX.xaxis.set_major_formatter(nullfmt)
axHistY.yaxis.set_major_formatter(nullfmt)

# Draw scatter plot
# , 'r*', mew=0, linewidth=0, linestyle='', marker='*', markerfacecolor='r', markersize=8, label='N(H) from Tau353'
axScatter.errorbar(x, y, xerr=xer, yerr=yer,  color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='$R=N_{H-tau}/N_{H-fromRad}$') 
axScatter.errorbar(x, yy, xerr=xer, yerr=yyer,  color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='$R=N_{H-EBV}/N_{H-fromRad}$') 
# axScatter.scatter(x, y/y,  marker='d', color = 'b', edgecolor='none', s=10, alpha=1)
# axScatter.scatter(x, yy, marker='s', color = 'k', edgecolor='none', s=10, alpha=1)
axScatter.legend(loc='upper right', fontsize=18)

# Draw y-axis histogram
axHistY.hist(y, num_y_bins, ec='red', fc='none', histtype='bar', orientation='horizontal')
axHistY.hist(yy, num_y_bins, ec='b', fc='none', histtype='bar', orientation='horizontal')

# Save figure
# fig.savefig('python_plot.pdf')

plt.grid(True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=15)
plt.tight_layout()
plt.show()



## Ratio vs CNM ##
x   = cnm
xer = cnmer

## For Plotting ##
min_y_data, max_y_data = np.min(y), np.max(y)
binsize                = 0.25
num_y_bins             = np.floor((max_y_data - min_y_data) / binsize)

# Axes definitions
nullfmt           = plt.NullFormatter()
left, width       = 0.1, 0.5
bottom, height    = 0.1, 0.8
bottom_h = left_h = left + width + 0.02

rect_scatter      = [left, bottom, width, height]
rect_histx        = [left, bottom_h, width, 0.4]
rect_histy        = [left_h, bottom, 0.3, height]

# Generate initial figure, scatter plot, and histogram quadrants
fig       = plt.figure(221, figsize=(16., 10.))

axScatter = fig.add_subplot(121, position=rect_scatter)
axScatter.set_xlabel('$N_{CNM} [10^{20} cm^{-2}]$', fontsize=35, fontweight='bold')
axScatter.set_ylabel('$R=N_{H-fromX}/N_{H-fromRad}$', fontsize=35, fontweight='bold')
# axScatter.set_xlim(0., 130.)
# axScatter.set_ylim(0.5, 4.)
axScatter.tick_params(axis='x', labelsize=20)
axScatter.tick_params(axis='y', labelsize=15)
axScatter.set_ylim(0.5, 4.)
axScatter.grid(True)

axHistY   = fig.add_subplot(122, position=rect_histy)
axHistY.set_xlim(0, 30.)
axHistY.set_ylim(0.5, 4.)
axHistY.set_xlabel('$Histogram$', fontsize=35, fontweight='bold')

# Remove labels from histogram edges touching scatter plot
# axHistX.xaxis.set_major_formatter(nullfmt)
axHistY.yaxis.set_major_formatter(nullfmt)

# Draw scatter plot
# , 'r*', mew=0, linewidth=0, linestyle='', marker='*', markerfacecolor='r', markersize=8, label='N(H) from Tau353'
axScatter.errorbar(x, y, xerr=xer, yerr=yer,  color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='$R=N_{H-tau}/N_{H-fromRad}$') 
axScatter.errorbar(x, yy, xerr=xer, yerr=yyer,  color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='$R=N_{H-EBV}/N_{H-fromRad}$') 
# axScatter.scatter(x, y/y,  marker='d', color = 'b', edgecolor='none', s=10, alpha=1)
# axScatter.scatter(x, yy, marker='s', color = 'k', edgecolor='none', s=10, alpha=1)
axScatter.legend(loc='upper right', fontsize=18)

# Draw y-axis histogram
axHistY.hist(y, num_y_bins, ec='red', fc='none', histtype='bar', orientation='horizontal')
axHistY.hist(yy, num_y_bins, ec='b', fc='none', histtype='bar', orientation='horizontal')

# Save figure
# fig.savefig('python_plot.pdf')

plt.grid(True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=15)
plt.show()