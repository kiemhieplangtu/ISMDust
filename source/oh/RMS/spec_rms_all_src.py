import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add fnewer of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import module            as md
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

## Velocity range ##
 #
 # return dict vel_info Infor for vel-range
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def vel_range_info():
	fname    = '../data/vel_range.txt'
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
def cal_rms(data, inf408, bd=2):
	### 93 MS sources ###
	new      = md.read_93src_info(fname = '../../hi/result/nhi_thin_cnm_wnm_93src.txt')
	vel_info = vel_range_info()
 	src_list = list(data.la.srcname)

 	print 'src;    l; b;   e(-t)_RMS;  tau_rms;    EM_RMS'
 	arms     = []
 	erms     = []
 	tau_sig  = {}
 	for src in src_list:
 		if(src not in new):
 			tau_sig[src] = -1
 			continue

		if(src == '3C223'): # no OH data for 3C223
			tau_sig[src] = -1
 			continue

	 	n,ell,bee,oh_f1,vlsr1,oh_f2,vlsr2,em_avg1,ab_avg1,em_avg2,ab_avg2 = get_src_info(data,src,src_list)
	 	avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2           = read_bins_to_cal_bg(n)
	 	xmin, xmax                                                        = vel_range(vel_info, n)

		# if (xmin == 0. and xmax == 0.):  # For OH-detected sources, Plot for Co-author paper
		# 	continue

	 	## Absorption and Emission data ##
	 	xd  = [vlsr1,           vlsr2]
		td  = [ab_avg1,         ab_avg2]
		tde = [em_avg1,         em_avg2]
		cst = [3.99757843817,   2.21841824609]
		frq = [1665.402,        1667.359]
		pfl = ['../data/gauss_1665_peaks.txt',   '../data/gauss_1667_peaks.txt']
			
		## OH1665 or OH1667 ##
		xd  = xd[bd-1]
		td  = td[bd-1]
		tde = tde[bd-1]
		cst = cst[bd-1]
		frq = frq[bd-1]
		pfl = pfl[bd-1]

		## rms of e(-tau) & Tau_SIG, and RMS of emission line
		tc1665,tc_er, tausig        = rms_baseline_from_linear_fit(xd, td, avmin1, avmax1, avmin2, avmax2,fit=False)
		bg_off1665,bgoff_er, nothig = rms_baseline_from_linear_fit(xd, tde, evmin1, evmax1, evmin2, evmax2,fit=False)
		tau_sig[src]                = tausig

		## rms of e(-tau) and RMS of emission line
		arms.append(tc_er)
		erms.append(bgoff_er)

		# etaud      = td/tc1665
		# taud       = -np.log(etaud)
		# tau_sigma  = get_tb_sigma(xd, taud, avmin1, avmax1, avmin2, avmax2)

 		print n, src, '\t', ell, '\t', bee, '\t', tc_er, '\t', tausig, '\t', bgoff_er  ## oh_spec_rms.txt

 	print ''
 	print 'RMS of e(-tau): ', arms
 	print 'Mean: ', np.mean(np.array(arms))
 	print ''
 	print 'RMS of emission line: ', erms
 	print 'Mean: ', np.mean(np.array(erms))
 	print ''

 	print 'RMS of Tau: ', tau_sig
 	print ''

 	for sc in new:
 		if(sc not in tau_sig):
 			tau_sig[sc] = -1
 		print tau_sig[sc]

#============== MAIN ==============#
data   = readsav('../data/makelines.sav')     # data.la
inf408 = readsav('../data/tb_408.sav')        # l_cntr, b_cntr, tb_408, Continuum at 408MHz
cal_rms(data, inf408, bd=1) # bd1=65, bd2=67

sys.exit()