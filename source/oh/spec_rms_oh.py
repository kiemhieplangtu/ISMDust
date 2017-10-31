import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import copy

from scipy.io.idl        import readsav
from restore             import restore
from mpfit               import mpfit
from scipy.integrate     import quad


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
def peak_info(source,fname='data/gauss_1665_peaks.txt'):
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
def read_tbg408_healpy(fname='result/bg408_to_compare.txt'):
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

## Read vel-range to calculate background ##
 #
 # params int n Order of the Source
 # params string fname File-name
 #
 # return dict ret Info of Vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def read_bins_to_cal_bg(n,fname='sub_data/bins_to_cal_bg.txt'):
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

	cond   = avmin1[n] == avmax1[n] == avmin2[n] == avmax2[n] == evmin1[n] == evmax1[n] == evmin2[n] == evmax2[n]
	if(cond):
		avmin1[n] = evmin1[n] = -40.
		avmax1[n] = evmax1[n] = -20.
		avmin2[n] = evmin2[n] = 10.
		avmax2[n] = evmax2[n] = 30. 

	return avmin1[n],avmax1[n],avmin2[n],avmax2[n],evmin1[n],evmax1[n],evmin2[n],evmax2[n]

## Read vel-range to calculate background and RMS ##
 #
 # params 1-D array xd x-data
 # params 1-D array td T-data
 # params float vmin1 Range of velocity
 # params float vmax1 Range of velocity
 # params float vmin2 Range of velocity
 # params float vmax2 Range of velocity
 #
 # return slope and intercept of the Linear fit
 #        RMS and tau_RMS (Absorption line only)
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def rms_baseline_from_linear_fit(xd,td,vmin1,vmax1,vmin2,vmax2,fit=True):
	tb = []
	v  = []
	count = 0
	for i in range(0,len(xd)):
		if (((xd[i] > vmin1) and (xd[i] < vmax1)) or ((xd[i] > vmin2) and (xd[i] < vmax2))):
			tb.append(td[i])
			v.append(xd[i])
			count = count + 1

	if(fit):
	 	slope, bsline = np.polyfit(v, tb, 1)
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

	sigma   = np.sqrt( sigma/(count-1) )
	tau_sig = np.std( -np.log( np.array(tb) ) )
	# print '    ', sigma, sigma, np.std(tb) ## sigma = np.std(tb) => Good!! they're std

 	return bsline, sigma, tau_sig

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
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
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

## Velocity range ##
 #
 # return dict vel_info Infor for vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def vel_range_info():
	fname    = 'data/vel_range.txt'
	cols     = ['idx','src','vmin','vmax']
	fmt      = ['i','s','f','f']
	vel      = restore(fname, 2, cols, fmt)
	vel_info = vel.read()

	return vel_info

## Velocity range ##
 #
 # params int n order of Source
 # params dict vel_info Infor for vel-range
 #
 # return float xmin, xmax Vel-range to fit
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def vel_range(vel_info, n):
	vmin = vel_info['vmin']
	vmax = vel_info['vmax']	
	
	xmin = vmin[n]
	xmax = vmax[n]

	return xmin, xmax


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

## Compute RMS of tau and Emission line for OH1665, OH16667 line ##
 #
 # params dict data     Data
 # params dict inf408   Info about the Tb_background at 408MHz
 # params int bd        OH65 or OH67
 #
 # return RMS of tau and emission line
 #
 # version 02/2017 
 # author Nguyen Van Hiep ##
def cal_rms(data, inf408, bd=1):
	bg408    = read_tbg408_healpy(fname='result/bg408_to_compare.txt')
	vel_info = vel_range_info()
 	src_list = list(data.la.srcname)

 	print 'src;    l; b;   e(-t)_RMS;  tau_rms;    EM_RMS'
 	arms = []
 	erms = []
 	for src in src_list:
	 	n,ell,bee,oh_f1,vlsr1,oh_f2,vlsr2,em_avg1,ab_avg1,em_avg2,ab_avg2 = get_src_info(data,src,src_list)
	 	avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2           = read_bins_to_cal_bg(n)
	 	xmin, xmax                                                        = vel_range(vel_info, n)

	 	# if (src != '3C109'):
			# continue
		if (xmin == 0. and xmax == 0.):  # For OH-detected sources, Plot for Co-author paper
			continue

	 	## Absorption and Emission data ##
	 	xd  = [vlsr1,           vlsr2]
		td  = [ab_avg1,         ab_avg2]
		tde = [em_avg1,         em_avg2]
		cst = [3.99757843817,   2.21841824609]
		frq = [1665.402,        1667.359]
		pfl = ['data/gauss_1665_peaks.txt',   'data/gauss_1667_peaks.txt']
			
		## OH1665 or OH1667 ##
		xd  = xd[bd-1]
		td  = td[bd-1]
		tde = tde[bd-1]
		cst = cst[bd-1]
		frq = frq[bd-1]
		pfl = pfl[bd-1]

		# xd,td = bin_up(xd,td,nbin=1)

		## Background ##
		tbg1665    = 2.8+get_tb_408(ell,bee,inf408.tb_408)*(408./frq)**2.8 # Tbg from 408MHz
		tbg1665    = 2.8+bg408[src]*(408./frq)**2.8 # Tbg from 408MHz
		if(src=='3C123'):
			tbg1665 = 25.0 #18.

		## BACKGROUNDS and THEIR UNCERTAINTIES ##
		## rms of e(-tau) & Tau_SIG, and RMS of emission line
		tc1665,tc_er, tau_sig       = rms_baseline_from_linear_fit(xd, td, avmin1, avmax1, avmin2, avmax2,fit=False)
		bg_off1665,bgoff_er, nothig = rms_baseline_from_linear_fit(xd, tde, evmin1, evmax1, evmin2, evmax2,fit=False)
		trx                         = bg_off1665 - tbg1665

		## rms of e(-tau) and RMS of emission line
		arms.append(tc_er)
		erms.append(bgoff_er)

		## 1-SIGMA STANDARD DEVIATION OF Tabspt and Temmission ##
		tab_sigma  = tc_er
		tem_sigma  = bgoff_er
		trx_sigma  = bgoff_er
		ton_sig    = np.sqrt(tem_sigma**2 + tab_sigma**2)

		## COMPUTE EXP(-TAU), TAU & 1-SIGMA STANDARD DEVIATION OF TAU ##
		etaud      = td/tc1665
		taud       = -np.log(etaud)
		tau_sigma  = get_tb_sigma(xd, taud, avmin1, avmax1, avmin2, avmax2)
		etau_sigma = np.abs(etaud)*np.sqrt( (tab_sigma/td)**2 + (tc_er/tc1665)**2 )
		etau_sigma = get_tb_sigma(xd, etaud, avmin1, avmax1, avmin2, avmax2)

 		print n, src, '\t', ell, '\t', bee, '\t', tc_er, '\t', tau_sig, '\t', bgoff_er  ## oh_spec_rms.txt

 	print ''
 	print 'RMS of e(-tau): ', arms
 	print 'Mean: ', np.mean(np.array(arms))
 	print ''
 	print 'RMS of emission line: ', erms
 	print 'Mean: ', np.mean(np.array(erms))
 	print ''

#============== MAIN ==============#
data   = readsav('data/makelines.sav')     # data.la
inf408 = readsav('data/tb_408.sav')        # l_cntr, b_cntr, tb_408, Continuum at 408MHz
cal_rms(data, inf408, bd=2) # bd1=65, bd2=67

sys.exit()