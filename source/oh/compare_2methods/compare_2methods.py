import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
import pylab             as pl
import oh_module         as ohmd
import module            as md

from restore             import restore


dat1    = ohmd.read_oh_carl(fname = '../../plotlib/oh/infor_oh_components.txt')
src11   = dat1['src']
tex11   = dat1['tex1']
tex12   = dat1['tex2']
tex11er = dat1['tex1er']
tex12er = dat1['tex2er']
noh12   = dat1['noh2']
noh12er = dat1['noh2er']
tau11   = dat1['tau1']
tau12   = dat1['tau2']
tau11er = dat1['tau1er']
tau12er = dat1['tau2er']
wid12   = dat1['wid2']
wdi12er = dat1['wid2er']

dat2    = ohmd.read_oh_hiep(fname = '../NOH/infor_oh_components_simple.txt')
src21   = dat2['src']
tex21   = dat2['tex1']
tex22   = dat2['tex2']
tex21er = dat2['tex1er']
tex22er = dat2['tex2er']
noh22   = dat2['noh2']
noh22er = dat2['noh2er']
tau21   = dat2['tau1']
tau22   = dat2['tau2']
tau21er = dat2['tau1er']
tau22er = dat2['tau2er']
wid22   = dat2['wid2']
wid22er = dat2['wid2er']

print len(tau12)

idx = 0
for i in range(len(src11)):
	if(tau12[i]==0 and tex12[i]==0):
		idx = i

tau12   = np.delete(tau12, idx)
tau12er = np.delete(tau12er,idx)
tex12   = np.delete(tex12,idx)
tex12er = np.delete(tex12er,idx)
tex22   = np.delete(tex22,idx)
tex22er = np.delete(tex22er,idx)
tau22   = np.delete(tau22,idx)
tau22er = np.delete(tau22er,idx)

print len(tau12)

# Plot
mpl.rcParams['axes.linewidth'] = 2.0
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

fig          = plt.figure(figsize=(14,10))
lbsize       = 22
ftsz         = 32

# Subplot 1
ax           = fig.add_subplot(121); #ax.set_rasterized(True)
major_xticks = np.arange(0.0, 0.4, 0.05)
minor_xticks = np.arange(0.0, 0.4, 0.025)
major_yticks = np.arange(0.0, 0.4, 0.05)
minor_yticks = np.arange(0.0, 0.4, 0.025)

ax.errorbar(tau11, tau21, xerr=tau11er, yerr=tau21er, color='gray', marker='o', ls='None', markersize=10, markeredgecolor='gray', markeredgewidth=1, label='data')
ax.errorbar(tau12, tau22, xerr=tau12er, yerr=tau22er, color='k', marker='o', ls='None', markersize=10, markeredgecolor='k', markeredgewidth=1, label='data')
ax.plot([0,0.3],[0,0.3], 'k--', label='')

ax.set_ylabel(r'$\mathrm{\tau_{_{PW}}}$', fontsize=ftsz+4,fontweight='normal')
ax.set_xlabel(r'$\mathrm{\tau_{_{Li(2017)}}}$', fontsize=ftsz+4, fontweight='normal')

ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
ax.tick_params(axis='x', labelsize=lbsize, pad=8)
ax.tick_params(axis='y', labelsize=lbsize, pad=2)
ax.tick_params(which='both', width=2)
ax.tick_params(which='major', length=12)
ax.tick_params(which='minor', length=8)

ax.set_ylim(0.0025, 0.3)
ax.set_xlim(0.0025, 0.3)
ax.grid(False)
ax.set_yscale('log')
ax.set_xscale('log')

## Subplot 2
ax1          = fig.add_subplot(122); #ax.set_rasterized(True)
major_xticks = np.arange(0.0, 35., 5.)
minor_xticks = np.arange(-3.0, 35., 1.)
major_yticks = np.arange(0.0, 35., 5.)
minor_yticks = np.arange(-3.0, 35., 1.)

ax1.errorbar(tex11, tex21, xerr=tex11er, yerr=tex21er, color='gray', marker='o', ls='None', markersize=10, markeredgecolor='gray', markeredgewidth=1, label='data')
ax1.errorbar(tex12, tex22, xerr=tex12er, yerr=tex22er, color='k', marker='o', ls='None', markersize=10, markeredgecolor='k', markeredgewidth=1, label='data')
ax1.plot([0,35.0],[0,35.0], 'k--', label='')

ax1.set_ylabel(r'$\mathrm{T_{ex_{PW}\ [K]}}$', fontsize=ftsz,fontweight='normal', labelpad=-12)
ax1.set_xlabel(r'$\mathrm{T_{ex_{Li(2017)}\ [K]}}$', fontsize=ftsz, fontweight='normal')

ax1.set_xticks(major_xticks)                                                       
ax1.set_xticks(minor_xticks, minor=True)                                           
ax1.set_yticks(major_yticks)                                                       
ax1.set_yticks(minor_yticks, minor=True)
ax1.tick_params(axis='x', labelsize=lbsize, pad=8)
ax1.tick_params(axis='y', labelsize=lbsize, pad=2)
ax1.tick_params(which='both', width=2)
ax1.tick_params(which='major', length=12)
ax1.tick_params(which='minor', length=8)
ax1.set_yscale('log')
ax1.set_xscale('log')

ax1.set_xlim(0.4, 30.0)
ax1.set_ylim(0.4, 30.0)
ax1.grid(False)

fig.subplots_adjust(wspace=20)

plt.tight_layout()
plt.savefig('tex_tau_compare.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=600)
plt.show()


sys.exit()






#### Histograms ###

fig = plt.figure(figsize=(8,6))
ax  = fig.add_subplot(111); ax.set_rasterized(True)
plt.subplot(1, 2, 1)
plt.hist(texratio, alpha=0.4, label='', color='k', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=20, lw=3)

plt.title('', fontsize=22)
plt.ylabel('Number of OH components', fontsize=22,fontweight='bold')
plt.xlabel('$Tex_{GD}/Tex_{CB}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=22)
plt.grid(False)


plt.subplot(1, 2, 2)
plt.hist(tauratio, alpha=0.4, label='', color='k', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=20, lw=3)

plt.title('', fontsize=22)
# plt.ylabel('Number of OH components', fontsize=22,fontweight='bold')
plt.xlabel(r'$\tau_{GD}/\tau_{CB}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=22)
plt.grid(False)

plt.tight_layout()
# plt.savefig('tex_tau_hist.eps', format='eps', dpi=600)
plt.show()

sys.exit()










## Tex Histogram ##
plt.hist(ratio, alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=40, lw=3)

plt.title('', fontsize=22)
plt.ylabel('Number of OH components', fontsize=22,fontweight='bold')
plt.xlabel('$Tex_{GD}/Tex_{CB}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=22)
plt.grid(False)
# plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(Hidx) (<3.0e20)', color='blue', fontsize=20)
plt.tight_layout()
plt.show()

## Tau67 Histogram ##
plt.hist(ratio, alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=40, lw=3)

plt.title('', fontsize=22)
plt.ylabel('Number of OH components', fontsize=22,fontweight='bold')
plt.xlabel(r'$\tau_{GD}/\tau_{CB}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=22)
plt.grid(False)
# plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(Hidx) (<3.0e20)', color='blue', fontsize=20)
plt.tight_layout()
plt.show()



## Tex67 - Linear line??? ##
plt.plot(tex12, tex22, label='', color='b', lw=3)
plt.plot([0,25], [0, 25], label='', color='k', lw=3)

plt.title('', fontsize=22)
plt.ylabel('Tex67_carl', fontsize=22,fontweight='bold')
plt.xlabel('Tex67_hiep', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=22)
plt.grid(False)
# plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(Hidx) (<3.0e20)', color='blue', fontsize=20)
plt.tight_layout()
plt.show()

## tau67 - Linear line??? ##
plt.plot(tau12, tau22, 'r.', ms=10)
plt.plot([0,0.25], [0, 0.25], label='', color='k', lw=3)

plt.title('', fontsize=22)
plt.ylabel('Tau67_carl', fontsize=22,fontweight='bold')
plt.xlabel('Tau67_hiep', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=22)
plt.grid(False)
# plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(Hidx) (<3.0e20)', color='blue', fontsize=20)
plt.tight_layout()
plt.show()

## Wid67 - Linear line??? ##
plt.plot(wid12, wid22, 'r.', ms=10)
plt.plot([0,2], [0, 2], label='', color='k', lw=3)

plt.title('', fontsize=22)
plt.ylabel('wid67_carl', fontsize=22,fontweight='bold')
plt.xlabel('wid67_hiep', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=22)
plt.grid(False)
# plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(Hidx) (<3.0e20)', color='blue', fontsize=20)
plt.tight_layout()
plt.show()



ret = {}
for i in range(len(src11)):
	if src11[i] not in ret.keys():
		ret[src11[i]]         = {}
		ret[src11[i]]['nohc'] = noh12[i]
		ret[src11[i]]['er2c'] = noh12er[i]**2

		ret[src11[i]]['nohh'] = noh22[i]
		ret[src11[i]]['er2h'] = noh22er[i]**2
	else:
		ret[src11[i]]['nohc'] = ret[src11[i]]['nohc'] + noh12[i]
		ret[src11[i]]['er2c'] = ret[src11[i]]['er2c'] + noh12er[i]**2

		ret[src11[i]]['nohh'] = ret[src11[i]]['nohh'] + noh22[i]
		ret[src11[i]]['er2h'] = ret[src11[i]]['er2h'] + noh22er[i]**2

		if( (i==len(src11)-1 ) or ( (i<len(src11)-1 ) and (src11[i] != src11[i+1])  ) ):
			ret[src11[i]]['er2h'] = np.sqrt( ret[src11[i]]['er2h'] )
			ret[src11[i]]['er2c'] = np.sqrt( ret[src11[i]]['er2c'] )

print ret['3C18']

noh_ratio = []
noh_c     = []
noh_h     = []
for sc in set(src11):
	print sc, ret[sc]['nohh']/ret[sc]['nohc']
	noh_ratio.append( ret[sc]['nohc']/ret[sc]['nohh'] )
	noh_c.append(ret[sc]['nohc'])
	noh_h.append(ret[sc]['nohh'])

## Histogram ##
plt.hist(noh_ratio, alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=40, lw=3)

plt.title('', fontsize=22)
plt.ylabel('Number of OH components', fontsize=22,fontweight='bold')
plt.xlabel('$NOH_{carl}/NOH_{hiep}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=22)
plt.grid(False)

# plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(Hidx) (<3.0e20)', color='blue', fontsize=20)

plt.tight_layout()
plt.show()

## Linear line??? ##
plt.plot(noh_c, noh_h, 'r.', label='')
plt.plot([0,8], [0, 8], label='', color='k', lw=3)

plt.title('', fontsize=22)
plt.ylabel('NOH67_carl', fontsize=22,fontweight='bold')
plt.xlabel('NOH67_hiep', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
# plt.legend(loc='upper right', fontsize=22)
plt.grid(False)

# plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(Hidx) (<3.0e20)', color='blue', fontsize=20)

plt.tight_layout()
plt.show()