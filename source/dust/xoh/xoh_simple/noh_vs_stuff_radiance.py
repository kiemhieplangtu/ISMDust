import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

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
def read_noh(fname = 'xoh_from_tau.txt'):
	cols = ['idx', 'src', 'l', 'b', 'xoh', 'xoher', 'nh2', 'nh2er', 'nhi', 'nhier', 'nh', 'nher', 'cnm', 'cnmer', 'noh', 'noher', 'av', 'aver']
	fmt  = ['i',   's',    'f','f', 'f'    , 'f'  , 'f'  , 'f'     , 'f'  , 'f'    , 'f'  , 'f'  , 'f'  , 'f'   , 'f'   , 'f'   , 'f'   , 'f' ]
	data = restore(fname, 4, cols, fmt)
	dat  = data.read(asarray=True)
	xoh  = dat['xoh']

	return dat['src'], dat['l'], dat['b'], dat['xoh'], dat['xoher'], dat['nh2'], dat['nh2er'], dat['nhi'], dat['nhier'], dat['nh'], dat['nher'], \
	dat['cnm'], dat['cnmer'], dat['noh'], dat['noher'], dat['av'], dat['aver']

#### MAIN #####
src, xl, xb, xoh, xoher, nh2, nh2er, nhi, nhier, nh, nher, cnm, cnmer, noh, noher, av, aver = read_noh(fname = 'noh_vs_stuff_radiance.txt')


## N(OH) vs Av ##
plt.errorbar(av, noh, xerr=aver, yerr=noher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.title('N$_{OH}$ (from Hiep) vs A$_{V}$', fontsize=30)
plt.xlabel('$A_{V}$ mag', fontsize=35)
plt.ylabel('$N_{OH} [10^{14} (cm^{-2}$]', fontsize=35)
plt.grid(True)
plt.ylim(-1., 8.)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)

# for i in range(len(src)):
# 	# if(oh[i] > 0):
# 	plt.annotate('('+str(src[i])+')', xy=(nh2[i], xoh[i]), xycoords='data',
#             xytext=(-50.,30.), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->"),fontsize=12,
#             )

plt.show()

## N(OH) vs NH2 ##
plt.errorbar(nh2, noh, xerr=nh2er, yerr=noher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.title('N$_{OH}$ (from Hiep) vs NH2', fontsize=30)
plt.xlabel('$NH2$', fontsize=35)
plt.ylabel('$N_{OH}$', fontsize=35)
plt.grid(True)
# plt.ylim(-10., 40.)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)

for i in range(len(src)):
	# if(oh[i] > 0):
	plt.annotate('('+str(src[i])+')', xy=(nh2[i], noh[i]), xycoords='data',
            xytext=(-50.,30.), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"),fontsize=12,
            )

plt.show()


## X(OH) vs NHI ##
plt.errorbar(nhi, noh, xerr=nhier, yerr=noher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.title('N$_{OH}$ (from Hiep) vs NHI', fontsize=30)
plt.xlabel('NHI', fontsize=35)
plt.ylabel('$N_{OH}$', fontsize=35)
plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)

for i in range(len(src)):
	# if(oh[i] > 0):
	plt.annotate('('+str(src[i])+')', xy=(nhi[i], noh[i]), xycoords='data',
            xytext=(-50.,30.), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"),fontsize=12,
            )
plt.show()

## X(OH) vs CNM ##
plt.errorbar(cnm, noh, xerr=cnmer, yerr=noher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.title('N$_{OH}$ (from Hiep) vs CNM', fontsize=30)
plt.xlabel('CNM', fontsize=35)
plt.ylabel('$N_{OH}$', fontsize=35)
plt.grid(True)
# plt.ylim(-10., 80.)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
for i in range(len(src)):
	# if(oh[i] > 0):
	plt.annotate('('+str(src[i])+')', xy=(cnm[i], noh[i]), xycoords='data',
            xytext=(-50.,30.), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"),fontsize=12,
            )
plt.show()

## N(OH) vs NH ##
plt.errorbar(nh, noh, xerr=nher, yerr=noher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.title('N$_{OH}$ (from Hiep) vs NH', fontsize=30)
plt.xlabel('NH', fontsize=35)
plt.ylabel('$X_{OH}$', fontsize=35)
plt.grid(True)
# plt.ylim(-10., 80.)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
for i in range(len(src)):
	# if(oh[i] > 0):
	plt.annotate('('+str(src[i])+')', xy=(nh[i], noh[i]), xycoords='data',
            xytext=(-50.,30.), textcoords='offset points',
            arrowprops=dict(arrowstyle="->"),fontsize=12,
            )
plt.show()