import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import module            as md

from numpy               import array
from restore             import restore

## Read info of XOH from tau, Ebv, Radiance #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 4/2017
 # Author Van Hiep ##
def read_nh_from_3proxies_for_OHsrc(fname = 'xoh_from_tau.txt'):
	cols = ['idx', 'src', 'l', 'b', 'xoh', 'xoher', 'nh2', 'nh2er', 'nhi', 'nhier', 'nh', 'nher', 'cnm', 'cnmer', 'noh', 'noher', 'av', 'aver']
	fmt  = ['i',   's',    'f','f', 'f'    , 'f'  , 'f'  , 'f'     , 'f'  , 'f'    , 'f'  , 'f'  , 'f'  , 'f'   , 'f'   , 'f'   , 'f'   , 'f' ]
	data = restore(fname, 4, cols, fmt)
	dat  = data.read(asarray=True)
	xoh  = dat['xoh']

	return dat['src'], dat['l'], dat['b'], dat['xoh'], dat['xoher'], dat['nh2'], dat['nh2er'], dat['nhi'], dat['nhier'], dat['nh'], dat['nher'], \
	dat['cnm'], dat['cnmer'], dat['noh'], dat['noher'], dat['av'], dat['aver']

#================= MAIN ========================#
t_src, t_xl, t_xb, t_xoh, t_xoher, t_nh2, t_nh2er, t_nhi, t_nhier, t_nh, t_nher, t_cnm, t_cnmer, t_noh, t_noher, t_av, t_aver = read_nh_from_3proxies_for_OHsrc(fname = 'xoh_from_tau.txt')
e_src, e_xl, e_xb, e_xoh, e_xoher, e_nh2, e_nh2er, e_nhi, e_nhier, e_nh, e_nher, e_cnm, e_cnmer, e_noh, e_noher, e_av, e_aver = read_nh_from_3proxies_for_OHsrc(fname = 'xoh_from_ebv2011.txt')
r_src, r_xl, r_xb, r_xoh, r_xoher, r_nh2, r_nh2er, r_nhi, r_nhier, r_nh, r_nher, r_cnm, r_cnmer, r_noh, r_noher, r_av, r_aver = read_nh_from_3proxies_for_OHsrc(fname = 'xoh_from_radiance.txt')

tyer = md.uncertainty_of_ratio(t_nh, r_nh, t_nher, r_nher)
eyer = md.uncertainty_of_ratio(e_nh, r_nh, e_nher, r_nher)
ryer = md.uncertainty_of_ratio(r_nh, r_nh, r_nher, r_nher)

### NH vs Av
plt.errorbar(t_av, t_nh/r_nh, xerr=t_aver, yerr=tyer, color='r', marker='o', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label='N(H) from Tau353')
plt.errorbar(e_av, e_nh/r_nh, xerr=e_aver, yerr=eyer, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='N(H) from E(B-V)')
plt.errorbar(r_av, r_nh/r_nh, xerr=r_aver, yerr=ryer, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='N(H) from Radiance')

plt.title('', fontsize=30)
plt.ylabel('$N_{H}\ [10^{20} cm^{-2}]$', fontsize=35, fontweight='bold')
plt.xlabel('$A_{v} [mag]$', fontsize=35, fontweight='bold')

plt.grid(True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc='upper left', fontsize=18)
plt.show()


sys.exit()

### NH vs NHI 
plt.errorbar(t_nhi, t_nh, xerr=t_nhier, yerr=t_nher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label='N(H) from Tau353')
plt.errorbar(e_nhi, e_nh, xerr=e_nhier, yerr=e_nher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='N(H) from E(B-V)')
plt.errorbar(r_nhi, r_nh, xerr=r_nhier, yerr=r_nher, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='N(H) from Radiance')

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


### NHI vs NH2 
plt.errorbar(t_nhi, t_nh2, xerr=t_nhier, yerr=t_nh2er, color='r', marker='o', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label='N(H) from Tau353')
plt.errorbar(e_nhi, e_nh2, xerr=e_nhier, yerr=e_nh2er, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='N(H) from E(B-V)')
plt.errorbar(r_nhi, r_nh2, xerr=r_nhier, yerr=r_nh2er, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='N(H) from Radiance')

plt.title('', fontsize=30)
plt.ylabel('$N_{H_{2}}\ [10^{20} cm^{-2}]$', fontsize=35, fontweight='bold')
plt.xlabel('$N_{HI} [10^{20} cm^{-2}]$', fontsize=35, fontweight='bold')

plt.grid(True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc='upper left', fontsize=18)
plt.show()


### NHI vs XOH
plt.errorbar(t_nhi, t_xoh, xerr=t_nhier, yerr=t_xoher, color='r', marker='o', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label='XOH from Tau353')
plt.errorbar(e_nhi, e_xoh, xerr=e_nhier, yerr=e_xoher, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='XOH from E(B-V)')
plt.errorbar(r_nhi, r_xoh, xerr=r_nhier, yerr=r_xoher, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='XOH from Radiance')

plt.title('', fontsize=30)
plt.ylabel('$X_{OH}$', fontsize=35, fontweight='bold')
plt.xlabel('$N_{HI} [10^{20} cm^{-2}]$', fontsize=35, fontweight='bold')

plt.grid(True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc='upper left', fontsize=18)
plt.show()


### XOH vs Av
plt.errorbar(t_av, t_xoh, yerr=t_xoher, xerr=t_aver, color='r', marker='o', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label='XOH from Tau353')
plt.errorbar(e_av, e_xoh, yerr=e_xoher, xerr=e_aver, color='b', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='XOH from E(B-V)')
plt.errorbar(r_av, r_xoh, yerr=r_xoher, xerr=r_aver, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='XOH from Radiance')

plt.title('', fontsize=30)
plt.ylabel('$X_{OH}$', fontsize=35, fontweight='bold')
plt.xlabel('$A_{V} [mag]$', fontsize=35, fontweight='bold')

plt.grid(True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc='upper left', fontsize=18)
plt.show()