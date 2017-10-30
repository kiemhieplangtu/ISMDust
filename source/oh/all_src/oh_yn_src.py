import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

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

	return tb_408[ibee, iell]

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

## Compute Tex for 1665 line ##
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

## Plot AVR & MED spectra ##
 #
 # params dict data     Data
 # params int  bd       OH665 or OH667
 #
 # return Tex and N(OH)
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def plot_avg_med(data, bd=1):
 	src_list = list(data.la.srcname)

 	for src in src_list:
	 	n       = src_list.index(src)

		ra50    = data.la.ra1950
		dec50   = data.la.dec1950
		ell     = data.la.ell
		bee     = data.la.bee

		oh_f1   = data.la.cfr_bd1
		vlsr1   = data.la.vlsr_bd1
		oh_f2   = data.la.cfr_bd2
		vlsr2   = data.la.vlsr_bd2

		em_avg1 = 0.5*correct_ctrl_chnl(data.la.i_em_avg_bd1)
		em_med1 = 0.5*correct_ctrl_chnl(data.la.i_em_med_bd1)
		ab_avg1 = 0.5*correct_ctrl_chnl(data.la.i_abs_avg_bd1)
		ab_med1 = 0.5*correct_ctrl_chnl(data.la.i_abs_med_bd1)

		em_avg2 = 0.5*correct_ctrl_chnl(data.la.i_em_avg_bd2)
		em_med2 = 0.5*correct_ctrl_chnl(data.la.i_em_med_bd2)
		ab_avg2 = 0.5*correct_ctrl_chnl(data.la.i_abs_avg_bd2)
		ab_med2 = 0.5*correct_ctrl_chnl(data.la.i_abs_med_bd2)

	 	avmin1,avmax1,\
	 	avmin2,avmax2,\
	 	evmin1,evmax1,\
	 	evmin2,evmax2 = read_bins_to_cal_bg(n)
	 	xmin, xmax    = vel_range(n)

	 	ohyn = 1
		if (xmin == 0. and xmax == 0.):
			ohyn = 0

	 	## Absorption and Emission data ##
	 	v    = vlsr1[n]
		t_ab = ab_avg1[n]
		t_em = em_avg1[n]
		y_ab = ab_med1[n]
		y_em = em_med1[n]
		cst  = 3.99757843817
		frq  = 1665.402
		pfl  = '../data/gauss_1665_peaks.txt'
		if(bd == 2):
			v    = vlsr2[n]
			t_ab = ab_avg2[n]
			t_em = em_avg2[n]
			y_ab = ab_med2[n]
			y_em = em_med2[n]
			cst  = 2.21841824609
			frq  = 1667.359
			pfl  = '../data/gauss_1667_peaks.txt'

		t_ab = t_ab/np.mean(t_ab)
		t_em = t_em/np.mean(t_em) + 0.02

		print ('{}    {}\t  {:08.4f}  {:08.4f}  {}'\
			.format(n, src, ell[n], bee[n], ohyn)) 


		continue

		## PLOT On-Source LINE ##
		plt.subplot(2,1,1)
		plt.plot(v,t_ab,'b-', label='AVG, On-source spectrum, $T_{on}$', ms=2)
		plt.plot(v,t_em,'r-', label='AVG, Off-source spectrum, $T_{off}$', ms=2)
		plt.title(src,fontsize=35)
		# plt.xlabel('$V_{lsr} (km/s)$',fontsize=35)
		plt.ylabel(r'$T (K)$',fontsize=35)
		plt.legend(loc='upper right')
		plt.grid()

		y_ab = y_ab/np.mean(y_ab)
		y_em = y_em/np.mean(y_em) + 0.02
		plt.subplot(2,1,2)
		plt.plot(v,y_ab,'b-', label='MED, On-source spectrum, $T_{on}$', ms=2)
		plt.plot(v,y_em,'r-', label='MED, Off-source spectrum, $T_{off}$', ms=2)
		# plt.title(src,fontsize=35)
		plt.xlabel('$V_{lsr} (km/s)$',fontsize=35)
		plt.ylabel(r'$T (K)$',fontsize=35)
		plt.legend(loc='upper right')
		plt.grid()

		plt.show()

#============== MAIN ==============#
data   = readsav('../data/makelines.sav') #data.la
plot_avg_med(data, bd=1)

sys.exit()