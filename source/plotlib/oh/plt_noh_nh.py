import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl

from restore             import restore

## Cal CNM uncertainties
 # Fit inf from Carl's paper
 #
 # params 
 #
 # return dict Infor for each src
 # 
 # version 2/2017
 # Author Van Hiep ##
def cal_cnm_uncertainty(tau, wid, tspin, sigtau, sigwid, sigtspin):
	d2  = wid**2   * tau**2    * sigtspin**2 + \
	      tspin**2 * wid**2    * sigtau**2 + \
	      tau**2   * tspin**2  * sigwid**2
	d   = 1.93988*np.sqrt(d2)

	return d

## Read infor for each HI src with associated OH components
 # Fit inf from Carl's paper
 #
 # inf void
 #
 # return dict Infor for each src
 # 
 # version 2/2017
 # Author Van Hiep ##
def read_nhi_with_assoc_oh_compnt(fname='data/asscociated_oh_hi_components.txt'):
	## Read fit inf of HI components from Carl papers ##
	cols = ['t_peak','err_t_peak','tau','err_tau','v0','err_v0','del_v','err_del_v', 'tspin', 'err_tspin', 'tkmax', 'nhi', 'cnm', 'frac', 'err_frac', 'src', 'oh_cpnt']
	fmt  = ['f',     'f',          'f',  'f',      'f', 'f',     'f',    'f',         'f',      'f',         'f',     'f',   'i',    'f',   'f',        's',  's'      ]
	dat  = restore(fname, 3, cols, fmt)
	inf  = dat.read()

	sc     = inf['src']
	tau    = inf['tau']
	tau_er = inf['err_tau']
	ts     = inf['tspin']
	ts_er  = inf['err_tspin']
	wid    = inf['del_v']
	wid_er = inf['err_del_v']
	cnm    = inf['cnm']

	nhi      = {}
	nhi_er   = {}
	oh_cpnt  = {}
	for i in range(len(sc)):
		if sc[i] not in nhi.keys():
			nhi[sc[i]]     = {}
			nhi_er[sc[i]]  = {}
			oh_cpnt[sc[i]] = {}

			xhi_er = 0.
			if(cnm[i] >= 0):
				xhi_er = cal_cnm_uncertainty(tau[i], wid[i], ts[i], tau_er[i], wid_er[i], ts_er[i])
				xhi_er = xhi_er/100.  ## xhi_er in [10^18] -> xhi_er in [10^20]

			nhi[sc[i]]     = [ inf['nhi'][i] ]
			nhi_er[sc[i]]  = [ xhi_er ]
			oh_cpnt[sc[i]] = [ inf['oh_cpnt'][i] ]
		else:
			xhi_er = 0.
			if(cnm[i] >= 0):
				xhi_er = cal_cnm_uncertainty(tau[i], wid[i], ts[i], tau_er[i], wid_er[i], ts_er[i])
				xhi_er = xhi_er/100. ## xhi_er in [10^18] -> xhi_er in [10^20]

			nhi[sc[i]]     += [ inf['nhi'][i] ]
			nhi_er[sc[i]]  += [ xhi_er ]
			oh_cpnt[sc[i]] += [ inf['oh_cpnt'][i] ]

	return nhi, nhi_er, oh_cpnt

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
	cols = ['id', 'src', 'l', 'b', 'tau1', 'tau1_er', 'v01', 'v01_er', 'wid1', 'wid1_er', 'tex1', 'tex1_er', 'noh1', 'noh1_er', 'ord1', 'tau2', 'tau2_er', 'v02', 'v02_er', 'wid2', 'wid2_er', 'tex2', 'tex2_er', 'noh2', 'noh2_er', 'ord2']
	fmt  = ['i',  's',  'f',  'f',  'f',   'f',        'f',  'f',      'f',    'f',       'f',     'f',      'f',     'f',      'i',     'f',    'f',       'f',   'f',      'f',    'f',       'f',     'f',      'f',    'f',       'i']
	dat  = restore(fname, 3, cols, fmt)
	inf  = dat.read()

	sc    = inf['src']
	noh1  = {}
	noh1e = {}

	noh2  = {}
	noh2e = {}
	for i in range(len(sc)):
		if sc[i] not in noh1.keys():
			noh1[sc[i]]  = {}
			noh1e[sc[i]] = {}

			noh2[sc[i]]  = {}
			noh2e[sc[i]] = {}

			noh1[sc[i]]  = [ inf['noh1'][i] ]
			noh1e[sc[i]] = [ inf['noh1_er'][i] ]

			noh2[sc[i]]  = [ inf['noh2'][i] ]
			noh2e[sc[i]] = [ inf['noh2_er'][i] ]
		else:
			noh1[sc[i]]  += [ inf['noh1'][i] ]
			noh1e[sc[i]] += [ inf['noh1_er'][i] ]

			noh2[sc[i]]  += [ inf['noh2'][i] ]
			noh2e[sc[i]] += [ inf['noh2_er'][i] ]

	return noh1, noh2, noh1e, noh2e

## Plot HI line and OH lines in same frame ##
 #
 # inf dict data Data
 #
 # return Void
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_hi_oh_assco_compnt():
	nhi, nhi_er, oh_cpnt       = read_nhi_with_assoc_oh_compnt(fname='data/associated_oh_hi_components.txt') ## Component-center of HI
	noh1, noh2, noh1e, noh2e   = read_noh_infor(fname='oh_components_infor.txt')  ## Component-center of OH1665, OH1667

	print "noh1, noh2, noh1e, noh2e"
	print noh1

	## e(-t)_rms, tau_rms, em_rms of OH to compute the Upper Limmit
	et_rms1, tau_rms1, em_rms1 = read_oh_rms_infor(fname='../../oh/result/oh1_spec_rms.txt') ## e(-t)_rms, tau_rms, em_rms of OH1665
	et_rms2, tau_rms2, em_rms2 = read_oh_rms_infor(fname='../../oh/result/oh2_spec_rms.txt') ## e(-t)_rms, tau_rms, em_rms of OH1667

	print tau_rms2

	## To compute the Upper Limmit
	tex1_mean                  = 3.  #6.02807 # K (See plt_tex_hist.py)
	tex2_mean                  = 3.5 #5.38855 # K
	wid_mean                   = 1.0 # km/s

	nhi_a   = []
	nhi_zr  = []
	nhi_ae  = []

	noh1_a  = []
	noh1_zr = []
	noh1_ae = []

	noh2_a  = []
	noh2_zr = []
	noh2_ae = []
	for sc in oh_cpnt:
		print sc
		oh_scpnt = oh_cpnt[sc]
		nhis     = nhi[sc]
		nhie_s   = nhi_er[sc]
		for i in range(len(nhis)):
			if( oh_scpnt[i] == '-' ):
				## Cal. upper limits from RMS, here use 3-sigma = 3*RMS
				nhi_zr.append(nhis[i])
				noh1_zr.append(4.26413*3.*tau_rms1[sc]*tex1_mean*wid_mean)
				noh2_zr.append(2.36634*3.*tau_rms2[sc]*tex2_mean*wid_mean)
			else:
				cpntlist = [int(x) for x in oh_scpnt[i].split(",")] ## One single HI component may have one/several OH components
				noh1_i   = 0.
				noh1_ie  = 0.

				noh2_i   = 0.
				noh2_ie  = 0.
				for j in range(len(cpntlist)):
					noh1_i  = noh1_i  + noh1[sc][cpntlist[j]]
					noh1_ie = noh1_ie + noh1e[sc][cpntlist[j]]**2

					noh2_i  = noh2_i  + noh2[sc][cpntlist[j]]
					noh2_ie = noh2_ie + noh2e[sc][cpntlist[j]]**2

				noh1_ie = np.sqrt(noh1_ie)
				noh2_ie = np.sqrt(noh2_ie)

				nhi_a.append(nhis[i])
				nhi_ae.append(nhie_s[i])
				noh1_a.append(noh1_i)
				noh1_ae.append(noh1_ie)
				noh2_a.append(noh2_i)
				noh2_ae.append(noh2_ie)

	# for i in range(len(nhi_a)):
	# 	print i, nhi_a[i], nhi_ae[i], noh1_a[i], noh2_a[i]

	## Plot
	mpl.rcParams['axes.linewidth'] = 1.5
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
	ax.errorbar(nhi_zr, noh2_zr, yerr=0.20*np.array(noh2_zr), color='gray', marker='v', ls='None', uplims=True, \
		markeredgewidth=1, markersize=9, capsize=2, markeredgecolor='gray', label='')
	xerb1, = plt.plot(nhi_zr, noh2_zr, color='gray', marker='v', ls='None', \
		markeredgewidth=1, markersize=9, markeredgecolor='gray', label=r'${\rm{Upper\ limit\ N}}_{\rm{OH}}(1667)$')

	ax.errorbar(nhi_a, noh2_a,xerr=nhi_ae, yerr=noh2_ae, \
		color='k', marker='o', ls='None', markersize=12, capsize=2, \
		markeredgecolor='k', markeredgewidth=1, label='')
	xerb2, = plt.plot(nhi_a, noh2_a, color='k', marker='o', ls='None', markersize=9, \
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
	leg   = plt.legend(fontsize=16, loc=(axbox.x0-0.09, axbox.y0+0.8))
	# leg.get_frame().set_linewidth(0.0)

	plt.grid(False)

	ax.set_xlim([0.015, 150.])
	ax.set_ylim([0.004, 9.])

	ax.set_xscale('log')
	ax.set_yscale('log')
	plt.tight_layout()	
	plt.savefig('Noh_vs_NHI.eps', format='eps', dpi=600)
	plt.show()

#============== MAIN ==============#
cal_hi_oh_assco_compnt()

sys.exit()