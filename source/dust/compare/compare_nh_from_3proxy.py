import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import module            as md

from numpy               import array
from restore             import restore

## Read NH obtained from 3 proxies for 35 los #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 08/2017
 # Author Van Hiep ##
def read_nh_from_3proxy(fname = 'nh_from_tau353.txt'):
	cols = ['idx','src','prx', 'nhi', 'nh','diff']
	fmt  = ['i',  's',  'f',   'f',   'f', 'f'   ]
	data = restore(fname, 3, cols, fmt)
	dat  = data.read(asarray=True)
	return dat

#================= MAIN ========================#
tauInfo = read_nh_from_3proxy(fname = 'nh_from_tau353.txt')
ebvInfo = read_nh_from_3proxy(fname = 'nh_from_ebv2011.txt')
radInfo = read_nh_from_3proxy(fname = 'nh_from_radiance.txt')

nhTau   = tauInfo['nh']
nhEbv   = ebvInfo['nh']
nhRad   = radInfo['nh']
indx    = radInfo['idx']
nhi     = radInfo['nhi']

plt.plot(nhi, nhTau, 'r*')
plt.plot(nhi, nhEbv, 'bh')
plt.plot(nhi, nhRad, 'kD')

plt.plot(nhi, nhTau, 'r*', mew=0, linewidth=0, linestyle='', marker='*', markerfacecolor='r', markersize=8, label='N(H) from Tau353')
plt.plot(nhi, nhEbv, 'bh', mew=0, linewidth=0, linestyle='', marker='h', markerfacecolor='b', markersize=8, label='N(H) from E(B-V)')
plt.plot(nhi, nhRad, 'kD', mew=0, linewidth=0, linestyle='', marker='D', markerfacecolor='k', markersize=8, label='N(H) from Radiance')

plt.title('', fontsize=30)
plt.ylabel('$N_{H}\ [10^{20} cm^{-2}]$', fontsize=35, fontweight='bold')
plt.xlabel('$N_{HI} [10^{20} cm^{-2}]$', fontsize=35, fontweight='bold')


plt.grid(True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc='upper left', fontsize=18)
plt.show()