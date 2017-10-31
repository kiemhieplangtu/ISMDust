import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import module            as md

from restore             import restore

## Read NHderived from 3 proxies from 35src HI only #
 #
 # params string fname Filename
 # return dict info 
 # 
 # version 8/2017
 # Author Van Hiep ##
def read_nh_from_proxies(fname = '../tau2nh/results/nh_from_ebv2011_35src_HI_only.dat'):
	cols = ['src', 'l', 'b', 'val', 'val_er', 'nh', 'nh_er', 'nhi', 'nhi_er', 'nh2', 'nh2_er']
	fmt  = ['s',   'f', 'f', 'f',   'f',      'f',  'f',     'f',   'f',      'f',   'f'   ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=True)

	return dat['val'], dat['val_er'], dat['nh'], dat['nh_er'], dat['nhi'], dat['nh_er']

####============= MAIN ==============###
Ebv, Ebver, Av, Src = md.read_ebv_av(fname = '../ebv2nh/data/ebv_sfd98_sf2011_for_35src_HI_only.txt', sfd98=False)
Aver                = 3.1*Ebver 

tau, tauer, tnh, tnher, nhi, nhier = read_nh_from_proxies(fname = '../tau2nh/results/nh_from_tau_35src_HI_only.dat')
ebv, ebver, enh, enher, nhi, nhier = read_nh_from_proxies(fname = '../ebv2nh/result/nh_from_ebv2011_35src_HI_only.dat')
r, rer, rnh, rnher, nhi, nhier     = read_nh_from_proxies(fname = '../radiance/results/nh_from_r_35src_HI_only.dat')

print nhi
print nhi.max()


# Plot - tau vs NH ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(10,5))
ax           = fig.add_subplot(111)

plt.errorbar(tau, tnh, xerr=tauer, yerr=tnher, color='b', marker='o', ls='-', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')
plt.errorbar(tau, nhi, xerr=tauer, yerr=nhier, color='k', marker='^', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')

plt.title('', fontsize=30)
plt.xlabel(r'$\mathrm{\tau_{353}}$', fontsize=22)
plt.ylabel(r'$\mathrm{N_{H} [10^{20} cm^{-2}]}$', fontsize=20)

plt.tick_params(axis='y', labelsize=18)
plt.tick_params(axis='x', labelsize=18)

plt.tick_params(axis='x', labelsize=16, pad=5)
plt.tick_params(axis='y', labelsize=16)
plt.tick_params(which='both', width=1.5)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)

plt.tight_layout()
# plt.savefig('qqq.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()









# Plot - R vs NH ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(10,5))
ax           = fig.add_subplot(111)

plt.errorbar(r, rnh, xerr=rer, yerr=rnher, color='b', marker='o', ls='-', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')
plt.errorbar(r, nhi, xerr=rer, yerr=nhier, color='k', marker='^', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')

plt.title('', fontsize=30)
plt.xlabel(r'$\mathrm{\mathcal{R}}$', fontsize=22)
plt.ylabel(r'$\mathrm{N_{H} [10^{20} cm^{-2}]}$', fontsize=20)

plt.tick_params(axis='y', labelsize=18)
plt.tick_params(axis='x', labelsize=18)

plt.tick_params(axis='x', labelsize=16, pad=5)
plt.tick_params(axis='y', labelsize=16)
plt.tick_params(which='both', width=1.5)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)

plt.tight_layout()
# plt.savefig('qqq.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()






# Plot - EBV vs NH ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(10,5))
ax           = fig.add_subplot(111)

plt.errorbar(ebv, enh, xerr=ebver, yerr=enher, color='b', marker='o', ls='-', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')
plt.errorbar(ebv, nhi, xerr=ebver, yerr=nhier, color='k', marker='^', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')

plt.title('', fontsize=30)
plt.xlabel(r'$\mathrm{E(B-V)}$', fontsize=22)
plt.ylabel(r'$\mathrm{N_{H} [10^{20} cm^{-2}]}$', fontsize=20)

plt.tick_params(axis='y', labelsize=18)
plt.tick_params(axis='x', labelsize=18)

plt.tick_params(axis='x', labelsize=16, pad=5)
plt.tick_params(axis='y', labelsize=16)
plt.tick_params(which='both', width=1.5)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)

plt.tight_layout()
# plt.savefig('qqq.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()











# Plot NH-ratios vs Av ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(10,5))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
# major_xticks = np.arange(1.e-7, 1.e-4, 1.e-5)
# minor_xticks = np.arange(1.e-7, 1.e-4, 5.e-6)
# major_yticks = np.arange(1., 80., 10.)                                              
# minor_yticks = np.arange(1., 80., 5.0)


plt.errorbar(Av, tnh/rnh, xerr=Aver, yerr=Aver*0., color='k', marker='o', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')
plt.errorbar(Av, enh/rnh, xerr=Aver, yerr=Aver*0., color='r', marker='^', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')
# plt.errorbar(nhi, rnh, xerr=nhier, yerr=rnher, color='b', marker='d', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')


plt.title('', fontsize=30)
plt.xlabel(r'$\mathrm{A_{V}}$', fontsize=22)
plt.ylabel(r'$\mathrm{N_{H}\ ratios}$', fontsize=20)

plt.tick_params(axis='y', labelsize=18)
plt.tick_params(axis='x', labelsize=18)
# plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))	

# ax.set_xticks(major_xticks)                                                       
# ax.set_xticks(minor_xticks, minor=True)                                           
# ax.set_yticks(major_yticks)                                                       
# ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=16, pad=5)
plt.tick_params(axis='y', labelsize=16)
plt.tick_params(which='both', width=1.5)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)

# plt.ylim(0.8, 90.0)
# plt.xlim(5.8e-7, 3.8e-5)
# plt.yscale('log')
# plt.xscale('log')

plt.tight_layout()
# plt.legend(loc='upper left', fontsize=18)
# plt.savefig('qqq.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()



# Plot NH-ratios vs NHI ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(10,5))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
# major_xticks = np.arange(1.e-7, 1.e-4, 1.e-5)
# minor_xticks = np.arange(1.e-7, 1.e-4, 5.e-6)
# major_yticks = np.arange(1., 80., 10.)                                              
# minor_yticks = np.arange(1., 80., 5.0)


plt.errorbar(nhi, tnh/rnh, xerr=nhier*0, yerr=nhier*0., color='k', marker='o', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='N(H) from tau')
plt.errorbar(nhi, enh/rnh, xerr=nhier*0, yerr=nhier*0., color='r', marker='^', ls='None', markersize=6, markeredgecolor='r', markeredgewidth=1, label='N(H) from E(B-V)')
# plt.errorbar(nhi, rnh, xerr=nhier, yerr=rnher, color='b', marker='d', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')
plt.axhline(1., xmin=0., xmax=40., ls='--', color='b', lw=2)


plt.title('', fontsize=30)
plt.xlabel(r'$\mathrm{N_{HI}}$', fontsize=22)
plt.ylabel(r'$\mathrm{N_{H}\ ratios}$', fontsize=20)

plt.tick_params(axis='y', labelsize=18)
plt.tick_params(axis='x', labelsize=18)
# plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))	

# ax.set_xticks(major_xticks)                                                       
# ax.set_xticks(minor_xticks, minor=True)                                           
# ax.set_yticks(major_yticks)                                                       
# ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=16, pad=5)
plt.tick_params(axis='y', labelsize=16)
plt.tick_params(which='both', width=1.5)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)

# plt.ylim(0.8, 90.0)
# plt.xlim(5.8e-7, 3.8e-5)
# plt.yscale('log')
# plt.xscale('log')

plt.tight_layout()
plt.legend(loc='upper left', fontsize=18)
# plt.savefig('qqq.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()

















# Plot NHI vs Av ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(10,5))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
# major_xticks = np.arange(1.e-7, 1.e-4, 1.e-5)
# minor_xticks = np.arange(1.e-7, 1.e-4, 5.e-6)
# major_yticks = np.arange(1., 80., 10.)                                              
# minor_yticks = np.arange(1., 80., 5.0)


plt.errorbar(nhi, Av, xerr=nhier, yerr=Av, color='k', marker='o', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')
# plt.errorbar(nhi, rnh, xerr=nhier, yerr=rnher, color='b', marker='d', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')


plt.title('', fontsize=30)
plt.ylabel(r'$\mathrm{A_{V}}$', fontsize=22)
plt.xlabel(r'$\mathrm{N_{HI} [10^{20} cm^{-2}]}$', fontsize=20)

plt.tick_params(axis='y', labelsize=18)
plt.tick_params(axis='x', labelsize=18)
# plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))	

# ax.set_xticks(major_xticks)                                                       
# ax.set_xticks(minor_xticks, minor=True)                                           
# ax.set_yticks(major_yticks)                                                       
# ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=16, pad=5)
plt.tick_params(axis='y', labelsize=16)
plt.tick_params(which='both', width=1.5)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)

# plt.ylim(0.8, 90.0)
# plt.xlim(5.8e-7, 3.8e-5)
# plt.yscale('log')
# plt.xscale('log')

plt.tight_layout()
# plt.legend(loc='upper left', fontsize=18)
# plt.savefig('qqq.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()


sys.exit()