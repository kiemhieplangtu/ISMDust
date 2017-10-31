import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import module            as md

from scipy.optimize      import curve_fit
from scipy.io.idl        import readsav
from numpy               import array
from restore             import restore
from plotting            import cplot

## Find the value of baseline ##
 #
 # inf list vlsr VLSR
 # inf list stock Stock-parameter

 # return float Value of baseline
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
def cal_bg(vlsr, stock):
	s     = 0.
	count = 0
	for i in range(0,2048):
		if (((vlsr[i] > -35.) and (vlsr[i] < -20.)) or ((vlsr[i] > 20.) and (vlsr[i] < 35.))):
			s     = s+stock[i]
			count = count + 1

	return s/count

## Interpolate the VLSR=0 bin ##
 #
 # inf list tb Temperature

 # return list tb Temperature
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
def correct_ctrl_chnl(tb):
	for i in range(0,101):
		tb[i][1023] = (tb[i][1021] + tb[i][1022] + tb[i][1024] + tb[i][1025])/4.

	return tb

## Get the index of a given velocity ##
 #
 # inf list v_axis Velocity axis
 # inf float vel Value of velocity
 #
 # return int idx Index of Velocity in the Vel_list
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

## Read infor for each HI src
 # Fit inf from Carl's paper
 #
 # inf void
 #
 # return dict Infor for each src
 # 
 # version 2/2017
 # Author Van Hiep ##
def read_nhi_v0_infor(fname='../hi/result/component_fit_params.txt'):
	## Read fit inf of HI components from Carl papers ##
	cols = ['t_peak','err_t_peak','tau','err_tau','v0','err_v0','del_v','err_del_v', 'tspin', 'err_tspin', 'tkmax', 'nhi', 'cnm', 'frac', 'err_frac', 'src']
	fmt  = ['f',     'f',          'f',  'f',      'f', 'f',     'f',    'f',         'f',      'f',         'f',     'f',   'i',    'f',   'f',        's']
	dat  = restore(fname, 3, cols, fmt)
	inf  = dat.read()

	sc   = inf['src']
	v0   = {}
	wid  = {}
	for i in range(len(sc)):
		if sc[i] not in v0.keys():
			v0[sc[i]]  = {}
			wid[sc[i]] = {}
			v0[sc[i]]  = [ inf['v0'][i] ]
			wid[sc[i]] = [ inf['del_v'][i] ]
		else:
			v0[sc[i]]  += [ inf['v0'][i] ]
			wid[sc[i]] += [ inf['del_v'][i] ]

	return v0,wid

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
	d2  = wid**2   * sigtspin**2 * tau + \
	      wid**2   * tspin**2    * sigtau**2 + \
	      tspin**2 * tau**2      * sigwid**2
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
def read_nhi_with_assoc_oh_compnt(fname='../plotlib/oh/data/asscociated_oh_hi_components.txt'):
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
			if(cnm[i] > 0):
				xhi_er = cal_cnm_uncertainty(tau[i], wid[i], ts[i], tau_er[i], wid_er[i], ts_er[i])

			nhi[sc[i]]     = [ inf['nhi'][i] ]
			nhi_er[sc[i]]  = [ xhi_er ]
			oh_cpnt[sc[i]] = [ inf['oh_cpnt'][i] ]
		else:
			xhi_er = 0.
			if(cnm[i] > 0):
				xhi_er = cal_cnm_uncertainty(tau[i], wid[i], ts[i], tau_er[i], wid_er[i], ts_er[i])

			nhi[sc[i]]     += [ inf['nhi'][i] ]
			nhi_er[sc[i]]  += [ xhi_er ]
			oh_cpnt[sc[i]] += [ inf['oh_cpnt'][i] ]

	return nhi, nhi_er, oh_cpnt

## Read infor for each OH src
 # Fit inf from Carl's paper
 #
 # inf void
 #
 # return dict Infor for each src
 # 
 # version 2/2017
 # Author Van Hiep ##
def read_oh_v0_infor(fname='../plotlib/oh/oh_components_infor.txt'):
	## Read fit inf of HI components from Carl papers ##
	cols = ['id', 'src', 'l', 'b', 'tau1', 'tau1_er', 'v01', 'v01_er', 'wid1', 'wid1_er', 'tex1', 'tex1_er', 'noh1', 'noh1_er', 'ord1', 'tau2', 'tau2_er', 'v02', 'v02_er', 'wid2', 'wid2_er', 'tex2', 'tex2_er', 'noh2', 'noh2_er', 'ord2']
	fmt  = ['i',  's',  'f',  'f',  'f',   'f',        'f',  'f',      'f',    'f',       'f',     'f',      'f',     'f',      'i',     'f',    'f',       'f',   'f',      'f',    'f',       'f',     'f',      'f',    'f',       'i']
	dat  = restore(fname, 3, cols, fmt)
	inf  = dat.read()

	sc   = inf['src']
	v01  = {}
	v02  = {}
	w01  = {}
	w02  = {}
	for i in range(len(sc)):
		if sc[i] not in v01.keys():
			v01[sc[i]] = {}
			v02[sc[i]] = {}
			w01[sc[i]] = {}
			w02[sc[i]] = {}
			v01[sc[i]] = [ inf['v01'][i] ]
			v02[sc[i]] = [ inf['v02'][i] ]
			w01[sc[i]] = [ inf['wid1'][i] ]
			w02[sc[i]] = [ inf['wid2'][i] ]
		else:
			v01[sc[i]] += [ inf['v01'][i] ]
			v02[sc[i]] += [ inf['v02'][i] ]
			w01[sc[i]] += [ inf['wid1'][i] ]
			w02[sc[i]] += [ inf['wid2'][i] ]

	return v01, v02, w01, w02


## Read infor for each OH src
 # Fit inf from Carl's paper
 #
 # inf void
 #
 # return dict Infor for each src
 # 
 # version 2/2017
 # Author Van Hiep ##
def read_noh_infor(fname='../plotlib/oh/oh_components_infor.txt'):
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

## Gaussian fit for 1667 line ##
 #
 # inf dict data Data
 #
 # return fit-results
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def hi_lines(data):
 	src_list = list(data.la.srcname)

	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee

	oh_f = data.la.cfr_bd0
	vlsr = data.la.vlsr_bd0

	em_avg = correct_ctrl_chnl(data.la.i_em_avg_bd0)
	em_med = correct_ctrl_chnl(data.la.i_em_med_bd0)
	ab_avg = correct_ctrl_chnl(data.la.i_abs_avg_bd0)
	ab_med = correct_ctrl_chnl(data.la.i_abs_med_bd0)

	# Peaks of considered sources #
	peak = {}
	peak['3C18']     = {'0':[-0.013,-8.04, 0.7]}
	peak['3C109']    = {'0':[-0.0037, 9.22, 0.52],'1':[-0.0056, 10.06, 0.63]}
	peak['3C409']    = {'0':[-0.0292, 15.35,0.62]}
	peak['3C132']    = {'0':[-0.0057,7.794,0.518]}
	peak['P0531+19'] = {'0':[-0.0011,1.88,0.9],'1':[-0.00095,5.047,1.0]}

	for src in src_list:
		n = src_list.index(src)

		# VLSR #
		x = vlsr[n]

		# On-/Off-source Tb and Background Continuum #
		t_on   = 0.5*ab_avg[n]
		t_off  = 0.5*em_avg[n]
		bg_on  = cal_bg(x,t_on)
		bg_off = cal_bg(x,t_off)

		# e(-tau) and tau #
		etau   = t_on/bg_on
		tau    = -np.log(etau)

		# Plot #
		fig    = cplot()
		trace1 = fig.lines(x,etau,label='exp(-tau)',
				prop=dict(color='r',
					      linewidth=1,
					      linestyle='solid',
					      marker='o',
					      markerfacecolor='b',
					      markersize=0
				))

		data   = [trace1]
		layout = dict(title  = 'HI 21 cm - '+str(n)+'th: '+src+'  exp(-tau)',
					  title_fontsize=30,
		              grid   = True,
		              legend = dict(loc='upper left', fontsize=18),
		              xaxis  = dict(label='vlsr (km/s)',tick_size=18,fontsize=35,xlim=[-60.,60.]),
		              yaxis  = dict(label='exp(-tau)',tick_size=18,fontsize=35),
		              text   = [dict(loc=[-0.5,0.4],text='',color='blue',fontsize=17),
		              			dict(loc=[-0.5,0.31],text='',color='red',fontsize=19)
		              		   ],
		             )

		fig.iplot(data,layout)


## Plot HI line and OH lines in same frame ##
 #
 # inf dict data Data
 #
 # return Void
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def plot_hi_oh(data):
	v0_hi, wid_hi      = read_nhi_v0_infor(fname='../hi/result/component_fit_params.txt') ## Component-center of HI
	v01, v02, w01, w02 = read_oh_v0_infor(fname='../plotlib/oh/oh_components_infor.txt')  ## Component-center of OH1665, OH1667

	src_list = list(data.la.srcname)

	src   = data.la.srcname
	ra50  = data.la.ra1950
	dec50 = data.la.dec1950
	ell   = data.la.ell
	bee   = data.la.bee

	hi_f0 = data.la.cfr_bd0
	oh_f1 = data.la.cfr_bd1
	oh_f2 = data.la.cfr_bd2

	vlsr0 = data.la.vlsr_bd0
	vlsr1 = data.la.vlsr_bd1
	vlsr2 = data.la.vlsr_bd2

	em_avg0 = data.la.i_em_avg_bd0
	em_med0 = data.la.i_em_med_bd0
	ab_avg0 = data.la.i_abs_avg_bd0
	ab_med0 = data.la.i_abs_med_bd0

	em_avg1 = data.la.i_em_avg_bd1
	em_med1 = data.la.i_em_med_bd1
	ab_avg1 = data.la.i_abs_avg_bd1
	ab_med1 = data.la.i_abs_med_bd1

	em_avg2 = data.la.i_em_avg_bd2
	em_med2 = data.la.i_em_med_bd2
	ab_avg2 = data.la.i_abs_avg_bd2
	ab_med2 = data.la.i_abs_med_bd2

	em_avg2 = correct_ctrl_chnl(em_avg2)
	em_med2 = correct_ctrl_chnl(em_med2)
	ab_avg2 = correct_ctrl_chnl(ab_avg2)
	ab_med2 = correct_ctrl_chnl(ab_med2)

	# Vrange infor #
	fname    = 'data/vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()
	vmin     = vel_info['vmin']
	vmax     = vel_info['vmax']

	# 26 src with no CO #
	fname    = 'data/26src_no_co.txt'
	cols     = ['idx','src']
	fmt      = ['i','s']

	src_no_co = restore(fname, 2, cols, fmt)
	s26info   = src_no_co.read()
	s26src    = s26info['src']

	for src in src_list:
		print src
		# if (src != 'P0531+19'):
		# 	continue

		n    = src_list.index(src)

		if((src in v01.keys()) and (src in v0_hi.keys()) ) :
			cen0 = v0_hi[src]
			cen1 = v01[src]
			cen2 = v02[src]
			wid0 = wid_hi[src]
			wid1 = w01[src]
			wid2 = w02[src]
		else:
			cen0 = []
			cen1 = []
			cen2 = []
			wid0 = []
			wid1 = []
			wid2 = []

		print '===='
		print src
		print cen0
		print wid0
		print cen1
		print cen2

		# VLSR #
		x0 = vlsr0[n]
		x1 = vlsr1[n]
		x2 = vlsr2[n]

		# On-/Off-source Tb and Background Continuum #
		t_on0   = ab_avg0[n]
		t_off0  = em_avg0[n]
		bg_on0  = cal_bg(x0,t_on0)
		bg_off0 = cal_bg(x0,t_off0)

		t_on1   = ab_avg1[n]
		t_off1  = em_avg1[n]
		bg_on1  = cal_bg(x1,t_on1)
		bg_off1 = cal_bg(x1,t_off1)

		t_on2   = ab_avg2[n]
		t_off2  = em_avg2[n]
		bg_on2  = cal_bg(x2,t_on2)
		bg_off2 = cal_bg(x2,t_off2)


		# e(-tau) and tau #
		etau0  = t_on0/bg_on0
		tau0   = -np.log(etau0)

		etau1  = t_on1/bg_on1
		tau1   = -np.log(etau1)
		rms1, v_fltr1, T_fltr1 = md.get_1sigma_spec(x1,etau1, -10., 16., vmin=-40., vmax=34.)

		etau2  = t_on2/bg_on2
		tau2   = -np.log(etau2)
		rms2, v_fltr2, T_fltr2 = md.get_1sigma_spec(x2,etau2, -10., 16., vmin=-40., vmax=34.)

		etau0  = (etau0-1.)*0.025+0.06
		etau1  = (etau1-1.)+0.03
		etau2  = (etau2-1.)

		# Plot #
		plt.figure(figsize=(10,10))
		plt.plot(x0,etau0, ls='-', color='b', label=r'$e^{\tau} - HI$',       lw=1)
		plt.plot(x1,etau1, ls='-', color='r', label=r'$e^{\tau} - OH(1665)$', lw=1)
		plt.axhline(y=-2.*rms1+0.03, lw=1)
		plt.plot(x2,etau2, ls='-', color='k', label=r'$e^{\tau} - OH(1667)$', lw=1)
		plt.axhline(y=-2.*rms2, lw=1)

		# for i in range(len(cen0)):
		# 	plt.axvline(cen0[i], color='b')
		# 	plt.plot( [cen0[i]-wid0[i]/2, cen0[i]+wid0[i]/2], [0.-i*0.005,0.-i*0.005], 'b-', lw=2)
		# 	plt.text( cen0[i]+0.2+wid0[i]/2, 0.-i*0.005, str(cen0[i]))

		# for i in range(len(cen1)):
		# 	plt.axvline(cen1[i], color='r')
		# 	plt.plot( [cen1[i]-wid1[i]/2, cen1[i]+wid1[i]/2], [0.+i*0.005,0.+i*0.005], 'r-', lw=2)
		# 	plt.text( cen1[i]+0.2+wid1[i]/2, 0.+i*0.005, str(cen1[i]))

		# for i in range(len(cen2)):
		# 	plt.axvline(cen2[i], color='k')
		# 	plt.plot( [cen2[i]-wid2[i]/2, cen2[i]+wid2[i]/2], [0.01+2*i*0.005,0.01+2*i*0.005], 'k-', lw=2)
		# 	plt.text( cen2[i]+0.2+wid2[i]/2, 0.01+2*i*0.005, str(cen2[i]))

		plt.title(src, fontsize=30)
		plt.ylabel('', fontsize=35)
		plt.xlabel('VLSR (km/s)', fontsize=35)
		plt.tick_params(axis='x', labelsize=20)
		plt.tick_params(axis='y', labelsize=20)
		plt.tick_params(which='both', width=2)
		plt.legend(loc='upper left', fontsize=18)
		plt.grid(True)
		plt.xlim(-60.,60.)
		# plt.ylim(0.,20.)
		plt.savefig('hi_oh_abs_spec_'+src+'.png', bbox_inches='tight' , format='png', dpi=300)
		# plt.show()


## Plot HI line and OH lines in same frame ##
 #
 # inf dict data Data
 #
 # return Void
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_hi_oh_assco_compnt():
	nhi, nhi_er, oh_cpnt     = read_nhi_with_assoc_oh_compnt(fname='../plotlib/oh/data/associated_oh_hi_components.txt') ## Component-center of HI
	noh1, noh2, noh1e, noh2e = read_noh_infor(fname='../plotlib/oh/oh_components_infor.txt')  ## Component-center of OH1665, OH1667

	nhi_a   = []
	nhi_ae  = []

	noh1_a  = []
	noh1_ae = []

	noh2_a  = []
	noh2_ae = []
	for sc in oh_cpnt:
		oh_scpnt = oh_cpnt[sc]
		nhis     = nhi[sc]
		nhie_s   = nhi[sc]
		for i in range(len(nhis)):
			if( oh_scpnt[i] == '-' ):
				noh1_i  = 0.
				noh1_ie = 0.

				noh2_i  = 0.
				noh2_ie = 0.
			else:
				cpntlist = [int(x) for x in oh_scpnt[i].split(",")] 
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
	# 	print nhi_a[i], noh1_a[i], noh2_a[i]

	## Plot
	plt.rc('font', weight='bold')
	plt.rc('text', usetex=True)
	plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

	fig          = plt.figure(figsize=(8,11.2))
	ax           = fig.add_subplot(111)                                 
	major_xticks = np.arange(0., 55., 10.)                                              
	minor_xticks = np.arange(0., 55, 1)
	major_yticks = np.arange(0., 5.6, 1)                                              
	minor_yticks = np.arange(0., 5.6, 0.1) 

	ax.plot(nhi_a, noh1_a, 'k*', label=r'${\rm{N}}_{\rm{OH}}(1665)$', markersize=14)
	ax.plot(nhi_a, noh2_a, 'r*', label=r'${\rm{N}}_{\rm{OH}}(1667)$', markersize=13)

	# ax.errorbar(nhi_a, noh1_a,xerr=nhi_ae, yerr=noh1_ae, color='k', marker='o', ls='None', markersize=14, markeredgecolor='b', markeredgewidth=1, label=r'${\rm{N}}_{\rm{OH}}(1665)$')
	# ax.errorbar(nhi_a, noh2_a,xerr=nhi_ae, yerr=noh2_ae, color='r', marker='o', ls='None', markersize=13, markeredgecolor='b', markeredgewidth=1, label=r'${\rm{N}}_{\rm{OH}}(1667)$')

	plt.title('', fontsize=22)
	plt.ylabel(r'${\rm{N}}_{\rm{OH}}$ $[10^{14} (cm^{-2})]$', fontsize=22)
	plt.xlabel(r'${\rm{N}}_{\rm{HI}}$ $[10^{20} (cm^{-2})]$', fontsize=22)
	ax.set_xticks(major_xticks)                                                       
	ax.set_xticks(minor_xticks, minor=True)                                           
	ax.set_yticks(major_yticks)                                                       
	ax.set_yticks(minor_yticks, minor=True)
	plt.tick_params(axis='x', labelsize=22, pad=8)
	plt.tick_params(axis='y', labelsize=22)
	plt.tick_params(which='both', width=2)
	plt.tick_params(which='major', length=9)
	plt.tick_params(which='minor', length=4)
	plt.legend(loc='upper right', fontsize=22)
	plt.grid(False)
	plt.xlim(0., 55.)
	plt.ylim(-0.2, 5.5)
	plt.tight_layout()
	# plt.savefig('Noh_vs_NHI.png', format='png', dpi=600)
	plt.show()

#============== MAIN ==============#
data = readsav('data/makelines.sav')
plot_hi_oh(data) #Plot HI & OH lines
# cal_hi_oh_assco_compnt()

sys.exit()