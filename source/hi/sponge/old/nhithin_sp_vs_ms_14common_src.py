import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import module            as md

from restore             import restore

## Read info of each source to dictionary #
 #
 # params dict dat Data
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read2dict(dat):
	sc     = dat['src']
	xl     = dat['l']
	xb     = dat['b']
	hi     = dat['nhi']
	hier   = dat['nhi_er']
	thin   = dat['thin']
	thiner = dat['thin_er']
	cnm    = dat['cnm']
	cnmer  = dat['cnm_er']
	wnm    = dat['wnm']
	wnmer  = dat['wnm_er']

	ret = {}
	for i in range(len(sc)):
		ret[sc[i]]           = {}
		ret[sc[i]]['l']      = xl[i]
		ret[sc[i]]['b']      = xb[i]
		ret[sc[i]]['nhi']    = hi[i]
		ret[sc[i]]['nhi_er'] = hier[i]
		ret[sc[i]]['thin']   = thin[i]
		ret[sc[i]]['thiner'] = thiner[i]
		ret[sc[i]]['cnm']    = cnm[i]
		ret[sc[i]]['cnm_er'] = cnmer[i]
		ret[sc[i]]['wnm']    = wnm[i]
		ret[sc[i]]['wnm_er'] = wnmer[i]

	return ret

## ======== MAIN ============= ##
## 21SPONGE Claire 30 Sources
sp30sc  = md.read_info_sponge_30src(fname = '../../oh/sponge/sub_data/30src_claire.txt')
sc1     = sp30sc['src']
xl1     = sp30sc['l']
xb1     = sp30sc['b']
hi1     = sp30sc['nhi']
hi1er   = sp30sc['nhi_er']
thin1   = sp30sc['thin']
thin1er = sp30sc['thin_er']
cnm1    = sp30sc['cnm']
cnm1er  = sp30sc['cnm_er']
wnm1    = sp30sc['wnm']
wnm1er  = sp30sc['wnm_er']
spdat   = read2dict(sp30sc)

## 78 MS Sources
ms78sc  = md.read_info_ms_78src(fname = '../result/nhi_lb_thin_cnm_wnm_78src.txt', scale=False)
sc2     = ms78sc['src']
xl2     = ms78sc['l']
xb2     = ms78sc['b']
hi2     = ms78sc['nhi']
hi2er   = ms78sc['nhi_er']
thin2   = ms78sc['thin']
thin2er = ms78sc['thin_er']
cnm2    = ms78sc['cnm']
cnm2er  = ms78sc['cnm_er']
wnm2    = ms78sc['wnm']
wnm2er  = ms78sc['wnm_er']
msdat   = read2dict(ms78sc)

## Common & different Sources
comsc = []  ## 14 src
difsc = []  ## 16 src
for sc in spdat:
	if (sc in msdat):
		comsc.append(sc)
	else:
		difsc.append(sc)

spnhi     = []
msnhi     = []
sp_er     = []
ms_er     = []
spwnm     = []
mswnm     = []
spwnmer   = []
mswnmer   = []
spcnm     = []
spcnmer   = []
mscnm     = []
mscnmer   = []
xl        = []
xb        = []
thinms    = []
thinms_er = []
thinsp    = []
thinsp_er = []
for i in range(len(comsc)):
	sc = comsc[i]
	# print sc, spdat[sc]['l'], msdat[sc]['l'], spdat[sc]['l'] - msdat[sc]['l']
	print sc, spdat[sc]['nhi'], msdat[sc]['nhi'], 100.*(spdat[sc]['nhi'] - msdat[sc]['nhi'])/msdat[sc]['nhi']
	spnhi.append(spdat[sc]['nhi'])
	msnhi.append(msdat[sc]['nhi'])
	spwnm.append(spdat[sc]['wnm'])
	mswnm.append(msdat[sc]['wnm'])
	spcnm.append(spdat[sc]['cnm'])
	spcnmer.append(spdat[sc]['cnm_er'])
	mscnm.append(msdat[sc]['cnm'])
	mscnmer.append(msdat[sc]['cnm_er'])

	spwnm.append(spdat[sc]['wnm'])
	spwnmer.append(spdat[sc]['wnm_er'])
	mswnm.append(msdat[sc]['wnm'])
	mswnmer.append(msdat[sc]['wnm_er'])

	ms_er.append(msdat[sc]['nhi_er'])
	sp_er.append(spdat[sc]['nhi_er'])
	xl.append(msdat[sc]['l'])
	xb.append(msdat[sc]['b'])
	thinms.append(msdat[sc]['thin'])
	thinms_er.append(msdat[sc]['thiner'])
	thinsp.append(spdat[sc]['thin'])
	thinsp_er.append(spdat[sc]['thiner'])

print 'NHI'
for i in range(len(comsc)):
	sc = comsc[i]
	print sc, spdat[sc]['nhi'], msdat[sc]['nhi'], 100.*(-spdat[sc]['nhi'] + msdat[sc]['nhi'])/msdat[sc]['nhi']

print 'WNM'
for i in range(len(comsc)):
	sc = comsc[i]
	print sc, spdat[sc]['wnm'], msdat[sc]['wnm'], 100.*(-spdat[sc]['wnm'] + msdat[sc]['wnm'])/msdat[sc]['wnm']

print 'CNM'
for i in range(len(comsc)):
	sc = comsc[i]
	if(msdat[sc]['cnm'] != 0.):
		print sc, spdat[sc]['cnm'], msdat[sc]['cnm'], 100.*(-spdat[sc]['cnm'] + msdat[sc]['cnm'])/msdat[sc]['cnm']
	else:
		print sc, spdat[sc]['cnm'], msdat[sc]['cnm']

########### MPFIT - N(HI) ############
spnhi = np.array(spnhi)
msnhi = np.array(msnhi)

spnhi_er = np.array(sp_er)
msnhi_er = np.array(ms_er)

xdata = np.array(msnhi)
ydata = np.array(spnhi)

# Error bar for x-axis and y-axis
xerr = np.array(msnhi_er)
yerr = np.array(spnhi_er)

xfit, yfit, mu, sig, m, ea = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=[1.0])

# Plot for paper ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(10,5))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
major_xticks = np.arange(0., 45., 10.)
minor_xticks = np.arange(0., 45., 5.)
major_yticks = np.arange(0., 55., 10.0)
minor_yticks = np.arange(0., 55., 5.0)

plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='k', marker='h', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='data')
plt.plot(xfit, mu, 'k-', mew=1, linewidth=2, marker='o', markerfacecolor='b', markersize=0, label='Linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
plt.plot([0,60],[0,60], 'k--', label='$N_{HI(SP)} = N_{HI(MS)}$')

plt.title('', fontsize=30)
plt.xlabel(r'$\mathrm{N_{HI} [10^{20} cm^{-2}] (MS)}$', fontsize=24, fontweight='normal')
plt.ylabel(r'$\mathrm{N_{HI} [10^{20} cm^{-2}] (SP)}$', fontsize=24)
# plt.ylim(-2.0, 55.0)
# plt.xlim(-2.0, 42.0)
plt.grid(False)
plt.tick_params(axis='y', labelsize=18)
plt.tick_params(axis='x', labelsize=15)
# plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
# plt.yscale('log')
# plt.xscale('log')

plt.text(0.0, 40.0, '$y = ['+str(m)+'\pm'+str(ea) +']\cdot x$', color='k', fontsize=20)

ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=16, pad=2)
plt.tick_params(axis='y', labelsize=16, pad=2)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.tight_layout()
# plt.legend(loc='upper left', fontsize=18)
# plt.savefig('nhi_scale_fct.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()
########### END - ODR ############




sys.exit()

















########### MPFIT - N(HI) thin ############
xdata = np.array(thinms)
ydata = np.array(thinsp)

# Error bar for x-axis and y-axis
xerr = np.array(thinms_er)
yerr = np.array(thinsp_er)

########### MPFIT fit ############
xfit, yfit, mu, sig, m, ea = md.do_linMPfit(xdata, ydata, xerr, yerr, lguess=[1.0])

# plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='$N^*_{HI-SP} vs N^{*}_{HI-MS}$')
# plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
# plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
# plt.plot([0,60],[0,60], 'k--', label='$N^*_{HI-SP} = N^{*}_{HI-MS}$')

# plt.title('$N^{*SP}_{HI}$ and $N^{*MS}_{HI}$ along 14 common lines-of-sight', fontsize=30)
# plt.ylabel('$N^{thin}_{HI} [10^{20}$ cm$^{-2}]$, from 21-SPONGE', fontsize=35)
# plt.xlabel('$N^{thin}_{HI} [10^{20}$ cm$^{-2}]$, MS', fontsize=35)
# plt.grid(True)
# plt.tick_params(axis='x', labelsize=18)
# plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

# plt.text(0.0, 40.0, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot (N^*_{HI}/10^{20})$', color='blue', fontsize=20)
# plt.legend(loc='upper left', fontsize=18)
# plt.show()
########### END - MPFIT ############




########### ODR fit ############
xdata = np.array(thinms)
ydata = np.array(thinsp)

# Error bar for x-axis and y-axis
xerr = np.array(thinms_er)
yerr = np.array(thinsp_er)

xfit, yfit, mu, sig, m, ea = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=[1.0])

# plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='$N^*_{HI(SP)}$ vs $N^{*}_{HI(MS)}$')
# plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
# plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
# plt.plot([0,60],[0,60], 'k--', label='$N^*_{HI(SP)} = N^{*}_{HI(MS)}$')

# plt.title('$N^{*SP}_{HI}$ and $N^{*MS}_{HI}$ along 14 common lines-of-sight', fontsize=30)
# plt.ylabel('$N^{*}_{HI} [10^{20}$ cm$^{-2}]$, from 21-SPONGE', fontsize=35)
# plt.xlabel('$N^{*}_{HI} [10^{20}$ cm$^{-2}]$, MS', fontsize=35)
# # plt.xlim(0.0, 2.0)
# # plt.ylim(-1.0, 6.0)
# plt.grid(True)
# plt.tick_params(axis='x', labelsize=18)
# plt.tick_params(axis='y', labelsize=15)
# # plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

# plt.text(0.0, 40.0, '$y = ['+str(m)+'\pm'+str(ea) +']\cdot x$', color='blue', fontsize=20)
# plt.legend(loc='upper left', fontsize=18)
# plt.show()

# Plot for paper ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(10,5))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
major_xticks = np.arange(0., 45., 10.)
minor_xticks = np.arange(0., 45., 5.)
major_yticks = np.arange(0., 55., 10.0)
minor_yticks = np.arange(0., 55., 5.0)

plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='k', marker='h', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, label='data')
plt.plot(xfit, mu, 'k-', mew=1, linewidth=2, marker='o', markerfacecolor='b', markersize=0, label='Linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
plt.plot([0,60],[0,60], 'k--', label='$N^*_{HI(SP)} = N^{*}_{HI(MS)}$')

plt.title('', fontsize=30)
plt.xlabel(r'$\mathrm{N^{*}_{HI} [10^{20} cm^{-2}] (MS)}$', fontsize=24, fontweight='normal')
plt.ylabel(r'$\mathrm{N^{*}_{HI} [10^{20} cm^{-2}] (SP)}$', fontsize=24)
plt.ylim(-2.0, 55.0)
plt.xlim(-2.0, 42.0)
plt.grid(False)
plt.tick_params(axis='y', labelsize=18)
plt.tick_params(axis='x', labelsize=15)
# plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
# plt.yscale('log')
# plt.xscale('log')

# plt.text(0.0, 40.0, '$y = ['+str(m)+'\pm'+str(ea) +']\cdot x$', color='k', fontsize=20)

ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=16, pad=2)
plt.tick_params(axis='y', labelsize=16, pad=2)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.tight_layout()
# plt.legend(loc='upper left', fontsize=18)
plt.savefig('nhi_scale_factor.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()
########### END - ODR ############




## Histogram N(HI) thin ##
plt.hist(ydata/xdata, alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=40, lw=3)

plt.title('Histogram of optically thin $N^{*}_{HI(SP)}$/$N^{*}_{HI(MS)}$ \n along 14 common lines-of-sight', fontsize=22)
plt.ylabel('Number', fontsize=22,fontweight='bold')
plt.xlabel('$N^{*}_{HI(SP)}/N^{*}_{HI(MS)}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=22)
plt.grid(False)

plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(HI) (<3.0e20)', color='blue', fontsize=20)

plt.tight_layout()
plt.show()