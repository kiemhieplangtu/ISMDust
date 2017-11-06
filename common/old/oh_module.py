import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import operator
import copy

from   restore           import restore
from   mpfit             import mpfit

## Fit the absorption line  ##
 #
 # params string src Source
 # params 1-D array xd x-data 
 # params 1-D array td y-data 
 # params 1-D array lguess Guess-parameters
 #
 # return Infos and Parameters of the fit
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def ab_fit(src,xd,td,lguess,xmin_id,xmax_id,evmin1,evmax1,evmin2,evmax2):
	npar    =  len(lguess)
	guessp  = np.array(lguess, dtype='float64')
	plimd   = [[False,False]]*npar
	plims   = [[0.,0.]]*npar
	parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]}
	pname   = ['bg']*npar
	pfix    = [False]*npar

	## Case: T0629+10 ##
	if(src == 'T0629+10'):
		pfix[2]  = True
		pfix[6]  = True
		# plimd[2] = [True,True]
		# plimd[6] = [True,True]
		# plims[2] = [guessp[3]-1.,guessp[3]+1.]
		# plims[6] = [guessp[6]-1.,guessp[6]+1.]

	count   = 1
	for i in range(npar):
		if(i==0):
			pname[i] = 'Tcont'

		if( ((i-1)%3 == 0) and (i>0) ):
			pname[i] = 'tau'+str(count)

		if( ((i-2)%3 == 0) and (i>0) ):
			pname[i] = 'v0'+str(count)

		if( ((i-3)%3 == 0) and (i>0) ):
			pname[i] = 'wid'+str(count)
			count    = count + 1

	parinfo = []
	for i in range(len(guessp)):
		parinfo.append(copy.deepcopy(parbase))

	for i in range(len(guessp)):
		parinfo[i]['value']   = guessp[i]
		parinfo[i]['fixed']   = pfix[i]
		parinfo[i]['parname'] = pname[i]
		parinfo[i]['limited'] = plimd[i]
		parinfo[i]['lims']    = plims[i]

	##  data and 1-sigma uncertainty ##
	x    = xd[xmin_id:xmax_id]
	y    = td[xmin_id:xmax_id]
	erry = np.zeros(y.shape, dtype=np.float64)
	sig  = get_tb_sigma(xd, td, evmin1, evmax1, evmin2, evmax2)
	erry.fill(sig)

	x    = x.astype(np.float64)
	y    = y.astype(np.float64)
	erry = erry.astype(np.float64)

	fa   = {'x':x, 'y':y, 'err':erry}
	mp   = mpfit(myabfunc, guessp, parinfo=parinfo, functkw=fa, quiet=True)

	## Fit values of Tau ##
	abp   = mp.params
	abper = mp.perror
	fit   = 0.
	for i in range(1, len(abp), 3):
		fit = fit + abp[i]*np.exp(- ( (x-abp[i+1])/(0.6005612*abp[i+2]))**2)  

	fit = np.exp(-fit)
	# fit = abp[1]*fit

	return x,fit,y/abp[0],abp,abper,npar,parbase,pname,parinfo

## Correct the central channel ##
 #
 # params 1D-array tb Temperature data
 # return dict ret All infor of line
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
def correct_ctrl_chnl(tb):
	for i in range(0,101):
		tb[i][1023] = (tb[i][1022] + tb[i][1024])/2.

	return tb

## Get the index of a given velocity ##
 #
 # params list v_axis Velocity axis
 # params float vel Value of velocity
 #
 # return int idx Index of Velocity in the Vel_list
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

## Read peaks Info ##
 #
 # params string fname File-name
 # return list guessp Guess-Parameters 
 #        list base_range Range of Baselines (in this script, don't care about it)
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def peak_info(source,fname='../data/gauss_1665_peaks.txt'):
	ret = {}

	cols = ['idx','src','tau','v0','wid','bmin','bmax']
	fmt  = ['i','s','f','f','f','f','f']
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read()
	bmin = info['bmin']
	bmax = info['bmax']
	src  = info['src']
	tau  = info['tau']
	v0   = info['v0']
	wid  = info['wid']

	for i in range(0,len(src)):
		if info['src'][i] not in ret.keys():
			ret[src[i]] = {}
			ret[src[i]]['base_range'] = [bmin[i],bmax[i]]
			ret[src[i]]['guessp']     = [tau[i],v0[i],wid[i]]
		else:
			ret[src[i]]['base_range'] = [bmin[i],bmax[i]]
			ret[src[i]]['guessp']     += [tau[i],v0[i],wid[i]]

	return ret[source]['guessp'], ret[source]['base_range']

## Create tau function ##
 #
 # params list x x axis
 # params list params Parameters
 #
 # return array y Multiple-Gaussian function
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def tau_func(x, *params):
	y = 0.
	for i in range(0, len(params), 4):
		amp = params[i]
		ctr = params[i+1]
		wid = params[i+2]
		tex = params[i+3]
		y   = y + tex*amp * np.exp( -((x - ctr)/(0.6005612*wid))**2)
	
	return y

## Retreive a SINGLE value of 408 t_b from haslam et al. ##
 #
 # params float ell Galactic-longitude
 # params float bee Galactic-latitude
 #
 # return float Tbg_408 Background-Temperature at (l,b)
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_tb_408(ell,bee,tb_408):
	iell= round(( modangle(ell)/360.)* 1080.)
	ibee= round( (( modangle( bee, 360., negpos=True)+90.)/180.)* 540)

	return tb_408[int(ibee), int(iell)]

## Convert angle to a specified range by finding the angle modulo the extent of the range. ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def modangle(angle, extent=360., negpos=False):
	offset = 0.
	if(negpos):
		offset = extent/2.

	return ((((angle-offset) % extent) + extent) % extent) - offset

## Read Tbg of 408MHz from Healpy map ##
 #
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def read_tbg408_healpy(fname='../result/bg408_to_compare.txt'):
	cols = ['src','l','b', 'il', 'ib', 'tbg', 'l-idx','b-idx','tbg1','tbg_hpy']
	fmt  = ['s','f','f','f','f','f','f','f','f','f']
	src  = restore(fname, 2, cols, fmt)
	info = src.read()

	src  = info['src']
	gl   = info['l']
	gb   = info['b']
	il   = info['il']
	ib   = info['ib']
	tbg  = info['tbg']
	lid  = info['l-idx']
	bid  = info['b-idx']
	bg1  = info['tbg1']
	bgh  = info['tbg_hpy']

	ret  = {}
	for i in range(len(src)):
		ret[src[i]] = bgh[i]

	return ret

## Read initial params to fit ##
 #
 # params string fname File-name
 #
 # return dict ret Info of Initial params
 #
 # version 07/2017 
 # author Nguyen Van Hiep ##
def read_init_params(fname='initial_info_to_fit.txt'):
	cols      = ['ord','src','vid1','vid2','tbg65','tbg67','tbgoff65','tbgoff67']
	fmt       = ['i',  's',  'i',    'i',  'f',     'f',  'f',         'f'      ]
	data      = restore(fname, 2, cols, fmt)
	dat       = data.read()
	
	src       = dat['src']
	vid1      = dat['vid1']
	vid2      = dat['vid2']

	tbg65     = dat['tbg65']
	tbg67     = dat['tbg67']

	tbgoff65  = dat['tbgoff65']
	tbgoff67  = dat['tbgoff67']

	ret       = {}
	for i in range(len(src)):
		# if src[i] not in ret.keys():
		ret[src[i]]             = {}
		ret[src[i]]['vid1']     = vid1[i]
		ret[src[i]]['vid2']     = vid2[i]

		ret[src[i]]['tbg65']    = tbg65[i]
		ret[src[i]]['tbg67']    = tbg67[i]

		ret[src[i]]['tbgoff65'] = tbgoff65[i]
		ret[src[i]]['tbgoff67'] = tbgoff67[i]

	return ret

## Read vel-range to calculate background ##
 #
 # params int n Order of the Source
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def read_bins_to_cal_bg(n,fname='../sub_data/bins_to_cal_bg.txt'):
	cols     = ['idx','src','avmin1','avmax1','avmin2','avmax2','evmin1','evmax1','evmin2','evmax2']
	fmt      = ['i','s','f','f','f','f','f','f','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()

	avmin1 = vel_info['avmin1']
	avmax1 = vel_info['avmax1']
	avmin2 = vel_info['avmin2']
	avmax2 = vel_info['avmax2']

	evmin1 = vel_info['evmin1']
	evmax1 = vel_info['evmax1']
	evmin2 = vel_info['evmin2']
	evmax2 = vel_info['evmax2']

	return avmin1[n],avmax1[n],avmin2[n],avmax2[n],evmin1[n],evmax1[n],evmin2[n],evmax2[n]

## Read vel-range to calculate background ##
 #
 # params 1-D array xd x-data
 # params 1-D array td T-data
 # params float vmin1 Range of velocity
 # params float vmax1 Range of velocity
 # params float vmin2 Range of velocity
 # params float vmax2 Range of velocity
 #
 # return slope and intercept of the Linear fit
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def baseline_from_linear_fit(xd,td,vmin1,vmax1,vmin2,vmax2,fit=True):
	if(fit):
		tb = []
		v  = []
		count = 0
		for i in range(0,len(xd)):
			if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
				tb.append(td[i])
				v.append(xd[i])
				count = count + 1

	 	slope,bsline = np.polyfit(v, tb, 1)
	else:
	 	bsline = 0.
		count  = 0.
		slope  = 0.
		for i in range(0,len(xd)):
			if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
				bsline = bsline + td[i]
				count  = count + 1

		bsline = bsline/count

 	sigma = 0.
 	for i in range(0,len(xd)):
		if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
			sigma = sigma + (td[i]-(slope*xd[i] + bsline))**2

	sigma = np.sqrt(sigma/(count-1))
 	return bsline, sigma

## Model of Absorption line and Emission line ##
 #
 # params 1-D array p Parameters
 # params fjac
 # params 1-D array x x-data
 # params 1-D array y y-data
 # params 1-D arrya err 1-sigma error of y-data
 #
 # return status of the fit
 #        weighted residuals
 #        Here use FWHM
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
# def tb_exp(x,y,bg,tau,v0,wid,tex,cont,err=1.):
def myabfunc(p, fjac=None, x=None, y=None, err=None):
	status  = 0
	ng      = (len(p)-1)/3

	fcn = 0.
	for i in range(1, len(p), 3):
		fcn = fcn + p[i]*np.exp(- ( (x-p[i+1])/(0.6005612*p[i+2]))**2)   

	exp_sumtau = np.exp(-fcn)
	tb_cont    = p[0]*exp_sumtau # P[0] is Tc

	return [status, (y-tb_cont)/err]

## Bin data-channel up by N neighbouring  points ##
 #
 # params x x-data
 # params y y-data
 # params n nbins
 #
 # return xdata, ydata
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def bin_up(x,t,nbin=4):
	chnl = 2048/nbin

	xx   = np.zeros(chnl, dtype=np.float64)
	yy   = np.zeros(chnl, dtype=np.float64)
	indx = 0
	for i in range(0,len(x),nbin):
		xtemp    = np.sum(x[i:(i+nbin)])
		ytemp    = np.sum(t[i:(i+nbin)])
		xx[indx] = xtemp/nbin
		yy[indx] = ytemp/nbin
		indx     =  indx + 1

	return xx,yy

## 1-sigma Error of tb-data ##
 #
 # params 1-D array xd x-data
 # params 1-D array td y-data
 # params float vmin1 vel-range
 # params float vmax1 vel-range
 # params float vmin2 vel-range
 # params float vmax2 vel-range
 #
 # return sigma
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_tb_sigma(xd, td, vmin1, vmax1, vmin2, vmax2):
	tb = []
	v  = []
	for i in range(len(xd)):
		if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
			tb.append(td[i])
			v.append(xd[i])

	v  = np.asarray(v, dtype=np.float64)
	tb = np.asarray(tb, dtype=np.float64)

 	return np.std(tb)

## Cal NOH uncertainties [10^14 cm-2]
 #
 # params float tau Value of tau
 # params float wid Value of wid
 # params float tspin Value of Tex
 # params float sigtau Error of tau
 # params float sigwid Error of FWHM-wid
 # params float sigtspin Error of Tex
 #
 # return Error of NOH1665
 # 
 # version 3/2017
 # Author Van Hiep ##
def cal_noh65_uncertainty(tau, wid, tspin, sigtau, sigwid, sigtspin):
	d2  = wid**2   * tau**2    * sigtspin**2 + \
	      tspin**2 * wid**2    * sigtau**2 + \
	      tau**2   * tspin**2  * sigwid**2
	return 4.26413*np.sqrt(d2) # [10^14 cm-2]

## Cal NOH uncertainties [10^14 cm-2]
 #
 # params float tau Value of tau
 # params float wid Value of wid
 # params float tspin Value of Tex
 # params float sigtau Error of tau
 # params float sigwid Error of FWHM-wid
 # params float sigtspin Error of Tex
 #
 # return Error of NOH1665
 # 
 # version 3/2017
 # Author Van Hiep ##
def cal_noh67_uncertainty(tau, wid, tspin, sigtau, sigwid, sigtspin):
	d2  = wid**2   * tau**2    * sigtspin**2 + \
	      tspin**2 * wid**2    * sigtau**2 + \
	      tau**2   * tspin**2  * sigwid**2
	return 2.36634*np.sqrt(d2) # [10^14 cm-2]

## Velocity range ##
 #
 # params int n order of Source
 #
 # return float xmin, xmax Vel-range to fit
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def vel_range(n):
	fname    = '../data/vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()
	vmin     = vel_info['vmin']
	vmax     = vel_info['vmax']	
	
	xmin = vmin[n]
	xmax = vmax[n]

	return xmin, xmax

# Read Tbg of 1666 from 408MHz #
#
# params string fname Filename
#
# return dict info of Tbg
# 
# Author Van Hiep
##
def read_tex_carl(fname = '../sub_data/tex_oh1165.txt'):
	cols = ['idx','src','amp','v0','wid','ts1','er1','ts2','er2','tbg', 'ts_carl', 'tau']
	fmt  = ['i','s','f','f','f','f','f','f','f','f','f','f']
	data = restore(fname, 3, cols, fmt)
	dat  = data.read()
	src  = dat['src']
	tau  = dat['tau']
	ts   = dat['ts_carl']

	ret  = {}
	for i in range(0,len(src)):
		if src[i] not in ret.keys():
			ret[src[i]] = {}
			ret[src[i]]['tau']     = [tau[i]]
			ret[src[i]]['ts_carl'] = [ts[i]]
		else:
			ret[src[i]]['tau']     += [tau[i]]
			ret[src[i]]['ts_carl'] += [ts[i]]

	return ret

## Get genaral info of the source ##
 #
 # params dict data Data
 # params str src Source-name
 #
 # return general Infos of source
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_src_info(data, src, src_list):
	n       = src_list.index(src)

	ra50    = data.la.ra1950
	dec50   = data.la.dec1950
	ell     = data.la.ell
	bee     = data.la.bee

	oh_f1   = data.la.cfr_bd1
	vlsr1   = data.la.vlsr_bd1
	oh_f2   = data.la.cfr_bd2
	vlsr2   = data.la.vlsr_bd2

	em_avg1 = correct_ctrl_chnl(data.la.i_em_avg_bd1)
	em_med1 = correct_ctrl_chnl(data.la.i_em_med_bd1)
	ab_avg1 = correct_ctrl_chnl(data.la.i_abs_avg_bd1)
	ab_med1 = correct_ctrl_chnl(data.la.i_abs_med_bd1)

	em_avg2 = correct_ctrl_chnl(data.la.i_em_avg_bd2)
	em_med2 = correct_ctrl_chnl(data.la.i_em_med_bd2)
	ab_avg2 = correct_ctrl_chnl(data.la.i_abs_avg_bd2)
	ab_med2 = correct_ctrl_chnl(data.la.i_abs_med_bd2)

	return n,ell[n],bee[n],oh_f1[n],vlsr1[n],oh_f2[n],vlsr2[n],em_avg1[n],ab_avg1[n],em_avg2[n],ab_avg2[n]

## Compute the Tex spectrum from etau ##
 #
 # params list etau_fit Exp(-tau)
 # params list velo     List of Velocity
 # params list tb_off   List of Off-source Temperature
 # params float tbg     CMB + Radio noise
 # params float trx     Receiver Temperature
 # params list popt     Fit-result parameters 
 #
 # return List tex      Excitation temperature
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def cal_tex_from_etau(etau_fit,velo,tb_off,tbg,trx,popt):
	tex  = []

	# Find the vel-range of the peak #
	peak = cal_vel_range_fwhm(popt)

	tex_peak  = []
	tex_count = []
	for k in range(0,len(peak)/2):
		tex_peak.append(0.)
		tex_count.append(0)

	# Compute the Excitation Temperature #
	s     = 0.
	count = 0
	for i in range(0, len(etau_fit)):
		if (etau_fit[i] != 1.):
			ts = (tb_off[i]-trx-tbg*etau_fit[i])/(1.-etau_fit[i])
			tex.append(ts)
		else:
			tex.append(0.)

		for k in range(0,len(peak),2):
			vmin = peak[0+k]
			vmax = peak[1+k]

			if ((velo[i]>vmin) and (velo[i]<vmax)) :
				tex_peak[k/2]  = tex_peak[k/2] + tex[i]
				tex_count[k/2] = tex_count[k/2] + 1

				s     = s + tex[i]
				count = count + 1

	print 'Tex for Peaks'
	ss = 0.
	for k in range(0,len(peak)/2):
		tex_peak[k]  = tex_peak[k]/tex_count[k]
		ss           = ss + tex_peak[k]
		print tex_peak[k],tex_count[k]
	print 'TEX ==== '
	print s/count, ss, count
	print 'TEX ==== '
	s = s/count
	return tex,s

## Compute the velocity-ranges within FHWM ##
 #
 # params list popt     Fit-result parameters 
 # return List intvl    Velocity-ranges within FHWM
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_peak_vel_range(popt):
	intvl = []
	for j in range(1, len(popt), 3):
		v0    = popt[1+j]                  # Line center
		wid   = popt[2+j]                  # Sigma of gaussian line
		intvl += [v0-wid/2., v0+wid/2.]    # vertical lines at 1/e

	return intvl


## Read oh-components from Carl ##
 #
 # params int n order of Source
 #
 # return infor
 #
 # version 07/2017 
 # author Nguyen Van Hiep ##
def read_oh_carl(fname = '../../plotlib/oh/infor_oh_components.txt'):
	cols = ['src','l','b','tau1','tau1er','v01','v01er','wid1','wid1er', 'tex1', 'tex1er', 'noh1', 'noh1er', 'tau2','tau2er','v02','v02er','wid2','wid2er', 'tex2', 'tex2er', 'noh2', 'noh2er']
	fmt  = ['s', 'f','f','f',    'f',    'f',   'f',   'f',   'f',      'f',    'f',       'f',     'f',     'f',    'f',    'f',   'f',   'f',   'f',      'f',    'f',       'f',     'f']
	data = restore(fname, 3, cols, fmt)
	dat  = data.read(asarray = True)
	return dat

## Read oh-components from Hiep ##
 #
 # params int n order of Source
 #
 # return infor
 #
 # version 07/2017 
 # author Nguyen Van Hiep ##
def read_oh_hiep(fname = '../NOH/infor_oh_components_simple.txt'):
	cols = ['src','l','b','tau1','tau1er','v01','v01er','wid1','wid1er', 'tex1', 'tex1er', 'noh1', 'noh1er', 'tau2','tau2er','v02','v02er','wid2','wid2er', 'tex2', 'tex2er', 'noh2', 'noh2er']
	fmt  = ['s', 'f','f','f',    'f',    'f',   'f',   'f',   'f',      'f',    'f',       'f',     'f',     'f',    'f',    'f',   'f',   'f',   'f',      'f',    'f',       'f',     'f']
	data = restore(fname, 3, cols, fmt)
	dat  = data.read(asarray = True)
	return dat