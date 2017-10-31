import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
import module            as md

from restore             import restore
from collections         import OrderedDict

## Read infor for each HI src with associated OH components
 # Fit inf from Carl's paper
 #
 # inf void
 #
 # return dict Infor for each src
 # 
 # version 2/2017
 # Author Van Hiep ##
def read_hi_oh_assoc_compnt(fname='data/asscociated_oh_hi_components.txt'):
	## Read fit inf of HI components from Carl papers ##
	cols  = ['t_peak','err_t_peak','tau','err_tau','v0','err_v0','del_v','err_del_v', 'tspin', 'err_tspin', 'tkmax', 'nhi', 'cnm', 'frac', 'err_frac', 'src', 'oh_cpnt']
	fmt   = ['f',     'f',          'f',  'f',      'f', 'f',     'f',    'f',         'f',      'f',         'f',     'f',   'i',    'f',   'f',        's',  's'      ]
	dat   = restore(fname, 3, cols, fmt)
	inf   = dat.read()

	nhi   = inf['nhi']
	tb    = inf['t_peak']
	tber  = inf['err_t_peak']

	tau   = inf['tau']
	tauer = inf['err_tau']

	v0    = inf['v0']
	v0er  = inf['err_v0']

	wid   = inf['del_v']
	wider = inf['err_del_v']

	Ts    = inf['tspin']
	Tser  = inf['err_tspin']

	tk    = inf['tkmax']
	cyn   = inf['cnm']

	src   = inf['src']

	ohcpt = inf['oh_cpnt']

	nhier = np.zeros(len(src))
	for i in range(len(src)):
		if(cyn[i] >= 0):
			nhier[i] = 0.01*md.get_NCNM_error(tau[i], tauer[i], Ts[i], Tser[i], wid[i], wider[i]) ## 1e18 -> 1e20
		else:
			nhier[i] = 0.01*md.get_NWNM_error(tb[i], tber[i], wid[i], wider[i])

	return src, nhi, nhier, ohcpt




## Read infor for each OH src
 # Fit inf from Carl's paper
 #
 # inf void
 #
 # return dict Infor for each src
 # 
 # version 2/2017
 # Author Van Hiep ##
def read_noh_infor(fname='oh_components_infor.txt'):
	## Read fit inf of HI components from Carl papers ##
	cols   = ['id', 'src', 'l', 'b', 'tau1', 'tau1_er', 'v01', 'v01_er', 'wid1', 'wid1_er', 'tex1', 'tex1_er', 'noh1', 'noh1_er', 'ord1', 'tau2', 'tau2_er', 'v02', 'v02_er', 'wid2', 'wid2_er', 'tex2', 'tex2_er', 'noh2', 'noh2_er', 'ord2']
	fmt    = ['i',  's',  'f',  'f',  'f',   'f',        'f',  'f',      'f',    'f',       'f',     'f',      'f',     'f',      'i',     'f',    'f',       'f',   'f',      'f',    'f',       'f',     'f',      'f',    'f',       'i']
	dat    = restore(fname, 3, cols, fmt)
	inf    = dat.read()

	idx    = inf['id']
	src    = inf['src']

	noh1   = inf['noh1']
	noh2   = inf['noh2']

	noh1er = inf['noh1_er']
	noh2er = inf['noh2_er']

	order  = inf['ord1']


	ret    = OrderedDict()
	for i in range(len(idx)):
		if(src[i] not in ret.keys()):
			ret[ src[i] ]           = {}
			ret[ src[i] ]['noh1']   = [ noh1[i] ]
			ret[ src[i] ]['noh2']   = [ noh2[i] ]
			ret[ src[i] ]['noh1er'] = [ noh1er[i] ]
			ret[ src[i] ]['noh2er'] = [ noh2er[i] ]
		else:
			ret[ src[i] ]['noh1']   += [ noh1[i] ]
			ret[ src[i] ]['noh2']   += [ noh2[i] ]
			ret[ src[i] ]['noh1er'] += [ noh1er[i] ]
			ret[ src[i] ]['noh2er'] += [ noh2er[i] ]

	return ret

## Cal N(OH) for associated HI component
 #
 # params list ohcpnt list of OH components associated with HI
 # params dict ohinfor Infor for OH
 #
 # return list [noh1, noh1er, noh2, noh2er]
 # 
 # version 2/2017
 # Author Van Hiep ##
def cal_assoc_noh(ohcpnt, ohinfor):
	noh1   = 0.
	noh1er = 0.
	noh2   = 0.
	noh2er = 0.
	for i in ohcpnt:
		noh1   = noh1   + ohinfor['noh1'][i]
		noh1er = noh1er + ohinfor['noh1er'][i]**2
		noh2   = noh2   + ohinfor['noh2'][i]
		noh2er = noh2er + ohinfor['noh2er'][i]**2

	noh1er = np.sqrt(noh1er)
	noh2er = np.sqrt(noh2er)

	return noh1, noh1er, noh2, noh2er


## Read RMS infor for each OH src
 #
 # params str fname Filename
 #
 # return dict RMS-Infor for each OH src
 # 
 # version 2/2017
 # Author Van Hiep ##
def read_oh_rms_infor(fname='../../oh/result/oh_spec_rms.txt'):
	## Read fit inf of HI components from Carl papers ##
	cols = ['id', 'src', 'l', 'b', 'et_rms', 'tau_rms', 'em_rms']
	fmt  = ['i',  's',   'f', 'f',  'f',      'f',        'f'   ]
	dat  = restore(fname, 2, cols, fmt)
	inf  = dat.read()

	sc       = inf['src']
	et_rms   = {}
	tau_rms  = {}
	em_rms   = {}
	for i in range(len(sc)):
			et_rms[ sc[i] ]  = inf['et_rms'][i]
			tau_rms[ sc[i] ] = inf['tau_rms'][i]
			em_rms[ sc[i] ]  = inf['em_rms'][i]

	return et_rms, tau_rms, em_rms


####--------------------MAIN----------------------###
## To compute the Upper Limmit
tex1_mean                  = 3.  #6.02807 # K (See plt_tex_hist.py)
tex2_mean                  = 3.5 #5.38855 # K
wid_mean                   = 1.0 # km/s


## e(-t)_rms, tau_rms, em_rms of OH to compute the Upper Limmit
et_rms1, tau_rms1, em_rms1 = read_oh_rms_infor(fname='../../oh/result/oh1_spec_rms.txt') ## e(-t)_rms, tau_rms, em_rms of OH1665
et_rms2, tau_rms2, em_rms2 = read_oh_rms_infor(fname='../../oh/result/oh2_spec_rms.txt') ## e(-t)_rms, tau_rms, em_rms of OH1667


src, \
nhi,\
nhier,\
ohcpnt = read_hi_oh_assoc_compnt(fname='data/associated_oh_hi_components.txt')

## All infor of OH 1665 and 1667
ohinfo = read_noh_infor(fname='oh_components_infor.txt')


print 'Src \t N(HI) \t N(OH67) '
noh1Cor   = np.zeros(len(src))
noh1erCor = np.zeros(len(src))
noh2Cor   = np.zeros(len(src))
noh2erCor = np.zeros(len(src))
nhiCor    = np.zeros(len(src))

nhiZro    = np.zeros(len(src))
noh1Zro   = np.zeros(len(src))
noh1erZro = np.zeros(len(src))
noh2Zro   = np.zeros(len(src))
noh2erZro = np.zeros(len(src))
for i in range(len(src)):
	sc    = src[i] 
	if(ohcpnt[i] != '-'):
		ohList    = [int(x) for x in ohcpnt[i].split(',')]
		nhiCor[i] = nhi[i]
		noh1Cor[i], \
		noh1erCor[i], \
		noh2Cor[i], \
		noh2erCor[i] = cal_assoc_noh(ohList, ohinfo[src[i]])

		print src[i], nhi[i], noh2Cor[i]
	else:
		nhiZro[i]    = nhi[i]
		noh1Zro[i]   = 4.26413*3.*tau_rms1[sc]*tex1_mean*wid_mean
		noh1erZro[i] = 0.
		noh2Zro[i]   = 2.36634*3.*tau_rms2[sc]*tex2_mean*wid_mean
		noh2erZro[i] = 0.


## Plot
mpl.rcParams['axes.linewidth'] = 2.0
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

fig          = plt.figure(figsize=(10,11))
ax           = fig.add_subplot(111)                                 
major_xticks = np.arange(0., 160., 20.)                                              
minor_xticks = np.arange(0., 160., 10.)
major_yticks = np.arange(-8., 9.0, 1)                                              
minor_yticks = np.arange(-8., 9.0, 0.1) 

# ax.plot(nhi_zr, noh1_zr, 'rv', label=r'${\rm{N}}_{\rm{OH}}(1665)$', markersize=8)
# ax.plot(nhi_zr, noh2_zr, 'kv', label=r'${\rm{N}}_{\rm{OH}}(1667)$', markersize=8)

## Plot Upper Limits
ax.errorbar(nhiZro, noh2Zro, yerr=0.20*np.array(noh2Zro), color='gray', marker='v', ls='None', uplims=True, \
markeredgewidth=1, markersize=9, capsize=2, markeredgecolor='gray', label='')
xerb1, = plt.plot(nhiZro, noh2Zro, color='gray', marker='v', ls='None', \
markeredgewidth=1, markersize=9, markeredgecolor='gray', label=r'${\rm{Upper\ limit\ N}}_{\rm{OH}}(1667)$')

ax.errorbar(nhi, noh2Cor,xerr=nhier, yerr=noh2erCor, \
color='k', marker='o', ls='None', markersize=12, capsize=2, \
markeredgecolor='k', markeredgewidth=1, label='')
xerb2, = plt.plot(nhi, noh2Cor, color='k', marker='o', ls='None', markersize=9, \
markeredgecolor='k', markeredgewidth=1, label=r'${\rm{N}}_{\rm{OH}}(1667)$')

plt.title('', fontsize=22)
plt.ylabel(r'${\rm{N}}_{\rm{OH}}$ $(10^{14} {\rm{cm}}^{-2})$', fontsize=22)
plt.xlabel(r'${\rm{N}}_{\rm{HI}}$ $(10^{20} {\rm{cm}}^{-2})$', fontsize=22)
ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)

axbox = ax.get_position()
leg   = plt.legend(fontsize=16, loc=(axbox.x0-0.099, axbox.y0+0.805))
# leg.get_frame().set_linewidth(0.0)

plt.grid(False)

ax.set_xlim([0.015, 150.])
ax.set_ylim([0.004, 9.])

ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()	
plt.savefig('Noh_vs_NHIa.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=150)
plt.show()