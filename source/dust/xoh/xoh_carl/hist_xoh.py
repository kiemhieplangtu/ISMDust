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
from scipy.odr           import *
## Read info of XOH from tau, Ebv, Radiance #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 4/2017
 # Author Van Hiep ##
def read_xoh(fname = 'xoh_from_tau.txt'):
	cols = ['idx','xoh']
	fmt  = ['i', 'f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	xoh  = dat['xoh']

	return np.array(xoh)

#### MAIN #####
t_xoh = read_xoh(fname = 'xoh_from_tau.txt')
e_xoh = read_xoh(fname = 'xoh_from_ebv2011.txt')
r_xoh = read_xoh(fname = 'xoh_from_radiance.txt')

print len(t_xoh)
print len(e_xoh)
print len(r_xoh)

c1           = 'r'
c2           = 'k' #'dimgrey'
c3           = 'mediumblue'
fig          = plt.figure(figsize=(8,11.2))
ax           = fig.add_subplot(111); ax.set_rasterized(True)

plt.hist(t_xoh, alpha=0.5,  label=r'${\rm{X_{OH},\ N_{H_{2}}\ from\ \tau_{353}}}$', color=c1, ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,90.0), bins=13, lw=3, edgecolor=c1)
plt.hist(e_xoh, alpha=0.5, label=r'${\rm{X_{OH},\ N_{H_{2}}\ from\ E(B-V)}}$', color=c2, ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,90.0), bins=13, lw=3, edgecolor=c2)
plt.hist(r_xoh, alpha=1.0,  label=r'${\rm{X_{OH},\ N_{H_{2}}\ from\ Radiance}}$',       color=c3, ls='--', histtype='step',       stacked=False, fill=False, range=(0.0,90.0), bins=13, lw=4, edgecolor=c3)

# major_ticks  = np.arange(0, 24., 5)
# minor_ticks  = np.arange(0, 24., 1)
# major_yticks = np.arange(0., 9., 2)
# minor_yticks = np.arange(0., 9., 1)
# plt.hist(t_xoh, alpha=0.5,  label=r'${\rm{X_{OH},\ N_{H_{2}}\ from\ \tau_{353}}}$', color=c1, ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,26.0), bins=13, lw=3, edgecolor=c1)
# plt.hist(e_xoh, alpha=0.5, label=r'${\rm{X_{OH},\ N_{H_{2}}\ from\ E(B-V)}}$', color=c2, ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,26.0), bins=13, lw=3, edgecolor=c2)
# plt.hist(r_xoh, alpha=1.0,  label=r'${\rm{X_{OH},\ N_{H_{2}}\ from\ Radiance}}$',       color=c3, ls='--', histtype='step',       stacked=False, fill=False, range=(0.0,26.0), bins=13, lw=4, edgecolor=c3)

# plt.hist(t_xoh, color=c1, ls='-',  histtype='stepfilled', stacked=False, fill=False, range=(0.0,26.0), bins=13, lw=3, edgecolor=c1)
# plt.hist(e_xoh, color=c2, ls='-',  histtype='stepfilled', stacked=False, fill=False, range=(0.0,26.0), bins=13, lw=3, edgecolor=c2)
# plt.hist(r_xoh, color=c3, ls='--', histtype='step',       stacked=False, fill=False, range=(0.0,26.0), bins=13, lw=4, edgecolor=c3)

# bins             = np.arange(0.,24.,2.)
# tex1, bin_edges1 = np.histogram(t_xoh, bins=bins, range=(0.0,24.), normed=False)
# tex2, bin_edges2 = np.histogram(e_xoh, bins=bins, range=(0.0,24.), normed=False)
# bins             = np.arange(0.,24.,2.)
# tbg,  bin_edges3 = np.histogram(r_xoh, bins=bins, range=(0.0,24.), normed=False)

# bincen1 = 0.5*(bin_edges1[1:] + bin_edges1[:-1])
# bincen2 = 0.5*(bin_edges2[1:] + bin_edges2[:-1])
# bincen3 = 0.5*(bin_edges3[1:] + bin_edges3[:-1])

# # plot1   = ax.plot(bincen1, tex1, ls='',  color=c1, marker='d',  markersize=12)
# # plot2   = ax.plot(bincen2, tex2, ls='',  color=c2, marker='d',  markersize=12)
# # plot3   = ax.plot(bincen3, tbg,  ls='',  color=c3, marker='d',  markersize=12)

# plot1   = ax.plot(bincen1, tex1, ls='',  color=c1, marker='*',  markersize=16)
# plot2   = ax.plot(bincen2, tex2, ls='',  color=c2, marker='d',  markersize=14)
# plot3   = ax.plot(bincen3, tbg,  ls='',  color=c3, marker='^',  markersize=14)

plt.title('OH abundance ratio: $X_{OH} = N_{OH}/N_{H_{2}}$', fontsize=22)
plt.ylabel('Number of sightlines', fontsize=22,fontweight='bold')
plt.xlabel('$X_{OH} = N_{OH}/N_{H_{2}} [10^{-8}]$', fontsize=22, fontweight='bold')
# ax.set_xticks(major_ticks)                                                       
# ax.set_xticks(minor_ticks, minor=True)                                           
# ax.set_yticks(major_yticks)                                                       
# ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=22)
plt.grid(False)
# plt.xlim(0.,23.)
# plt.ylim(0.,9.)

# plt.tight_layout()

# plt.savefig('xoh_hist.png', format='png', dpi=600)
plt.show()