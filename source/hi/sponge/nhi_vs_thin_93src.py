import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import module            as md
import copy

from   restore           import restore

## ======== MAIN ============= ##

## 93 MS+SPONGE Sources, 21SPONGE* prior
sc93   = md.read_nhi_93src(fname = '../result/nhi_thin_cnm_wnm_93src.txt')
sc     = []
xl     = []
xb     = []
hi     = []
hier   = []
thin   = []
thiner = []
cnm    = []
cnmer  = []
wnm    = []
wnmer  = []
for i in range(len(sc93['src'])):
	if(np.abs(sc93['b'][i]) > 0.):  ## exclude the Plane ???
		sc.append(sc93['src'][i])
		xl.append(sc93['l'][i])
		xb.append(sc93['b'][i])
		hi.append(sc93['nhi'][i])
		hier.append(sc93['nhi_er'][i])
		thin.append(sc93['thin'][i])
		thiner.append(sc93['thin_er'][i])

## 94 MS+SPONGE Sources, 21SPONGE* prior, R vs log10(N*HI)
fact   = []
lognhi = []
for i in range(0, len(thin)):
	temp = round(hi[i]/thin[i], 3)
	lognhi.append(np.log10(thin[i]*1.26))
	fact.append(temp)

# Error bar for x-axis and y-axis
thin   = np.asarray(thin)
hi     = np.asarray(hi)
xdata  = np.asarray(lognhi)
ydata  = np.asarray(fact)
hier   = np.array(hier)
thiner = np.array(thiner)
xerr   = thiner/thin/np.log(10.0)
yerr   = md.uncertainty_of_factors(ydata, thin, thiner, hi, hier)

print hi.min()
print hi.max()

########### MPFIT ############
xdata = xdata
ydata = ydata

print 'Length...'
print len(sc)

# Error bar for x-axis and y-axis
xerr = xerr
yerr = yerr

## Fit ##
# xfit, yfit, mu, sig, m, ea, b, eb = md.do_linMPfit(xdata, ydata, xerr, yerr, lguess=[ 0.33, 0.85])
xfit, yfit, mu, sig, m, ea, b, eb = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=[ 0.33, 0.85])  ## Both are fine and consitent

# Plot for paper ##
fig          = plt.figure(figsize=(10,5))
ax           = fig.add_subplot(111); #ax.set_rasterized(True) 

major_xticks = np.arange(0.0, 1.95, 0.2)
minor_xticks = np.arange(0.0, 1.95, 0.1)
major_yticks = np.arange(0.5, 2.25, 0.25)
minor_yticks = np.arange(0.5, 2.25, 0.125)

plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', mfc="None", markeredgewidth=1, capsize=0, label='data')
plt.plot(xfit, mu, 'k-', mew=1, linewidth=2, marker='o', markerfacecolor='b', markersize=0, label='Linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
# plt.plot([0,60],[0,60], 'k:', label='$N^*_{HI(SP)} = N^{*}_{HI(MS)}$')

ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=16, pad=6)
plt.tick_params(axis='y', labelsize=16, pad=2)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=6)
plt.tick_params(which='minor', length=2)

plt.title('', fontsize=30)
plt.ylabel('$f = N_{HI}$/$N^*_{HI}$', fontsize=22)
plt.xlabel('$log_{10}(N^*_{HI}/10^{20} cm^{-2}$)', fontsize=22)
plt.xlim(-0.0, 1.9)
plt.ylim(0.5, 2.0)
plt.grid(False)
plt.tick_params(axis='y', labelsize=18)
plt.tick_params(axis='x', labelsize=15)
# plt.yscale('log')
# plt.xscale('log')

# plt.text(0.1, 2.25, '$R = ['+str(m)+'\pm'+str(ea) +']\cdot log_{10}(N^*_{HI}/10^{20}) + ['+str(b)+'\pm'+str(eb)+']$', color='k', fontsize=18)
# plt.text(0.1, 2.15, r'$R = [0.32\pm0.06]\cdot log_{10}(N^*_{HI}/10^{20}) + [0.81\pm0.05]$, Lee et al. 2015', color='k', fontsize=18)

plt.tight_layout()
# plt.legend(loc='upper left', fontsize=18)
plt.show()


## For Plotting ##
min_y_data, max_y_data = np.min(ydata), np.max(ydata)
binsize                = 0.25
num_y_bins             = np.floor((max_y_data - min_y_data) / binsize)
num_y_bins             = 24

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
fig = plt.figure(1, figsize=(20, 10))
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

axScatter = plt.axes(rect_scatter)
# axHistx = plt.axes(rect_histx)
axHistY = plt.axes(rect_histy)

# Remove labels from histogram edges touching scatter plot
# axHistX.xaxis.set_major_formatter(nullfmt)
axHistY.yaxis.set_major_formatter(nullfmt)

major_xticks = np.arange(0.0, 2.0, 0.2)
minor_xticks = np.arange(0.0, 2.0, 0.1)
major_yticks = np.arange(0.75, 2.25, 0.25)
minor_yticks = np.arange(0.75, 2.25, 0.125)

# Draw scatter plot
flter = ydata<2.
axScatter.errorbar(xdata[flter], ydata[flter], xerr=xerr[flter], yerr=yerr[flter], color='k', marker='o', ls='None', markersize=5, markeredgecolor='k', markeredgewidth=1, capsize=0, label='data')
axScatter.plot(xfit, mu, 'k-', mew=1, linewidth=3, marker='o', markerfacecolor='b', markersize=0, label='Linear fit')
axScatter.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

axScatter.set_xticks(major_xticks)                                                       
axScatter.set_xticks(minor_xticks, minor=True)
axScatter.set_yticks(major_yticks)
axScatter.set_yticks(minor_yticks, minor=True)
axScatter.tick_params(axis='x', labelsize=26, pad=6)
axScatter.tick_params(axis='y', labelsize=24)
axScatter.tick_params(which='both', width=2)
axScatter.tick_params(which='major', length=14)
axScatter.tick_params(which='minor', length=8)

axScatter.set_ylim(0.75, 2.0)
axScatter.set_xlim(0.1, 1.9)
axScatter.grid(False)
# axScatter.legend(loc='upper left', fontsize=18)

# axHistY   = fig.add_subplot(122, position=rect_histy)
major_xticks = np.arange(0., 40., 10.)
minor_xticks = np.arange(0., 40., 2.)
major_yticks = np.arange(0., 3., 0.25)
minor_yticks = np.arange(0., 3., 0.125)

axScatter.set_xlabel(r'$\mathrm{log_{10}(N^*_{HI}/10^{20} cm^{-2}})$', fontsize=36, labelpad=13)
axScatter.set_ylabel(r'$f{=}\mathrm{N_{HI}/N^*_{HI}}$', fontsize=36)

# Draw y-axis histogram
axHistY.hist(ydata, num_y_bins, color='grey', ls='-',  histtype='stepfilled', stacked=False, fill=True, lw=1, orientation='horizontal')

plt.grid(False)
axHistY.set_xlabel(r'$\mathrm{\#\ of\ sightlines}$', fontsize=36, fontweight='normal', labelpad=5)

axHistY.set_xticks(major_xticks)                                                       
axHistY.set_xticks(minor_xticks, minor=True)
axHistY.set_yticks(major_yticks)
axHistY.set_yticks(minor_yticks, minor=True)
axHistY.tick_params(axis='x', labelsize=26, pad=6)
axHistY.tick_params(axis='y', labelsize=24)
axHistY.tick_params(which='both', width=2)
axHistY.tick_params(which='major', length=14)
axHistY.tick_params(which='minor', length=8)

axHistY.set_xlim(0, 38.)
axHistY.set_ylim(0.75, 2.)

# plt.tight_layout()
plt.savefig('ratio_vs_nhithin.eps', bbox_inches='tight', pad_inches=0.1, format='eps', dpi=600)
# plt.show()