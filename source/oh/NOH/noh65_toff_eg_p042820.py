import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl
import pylab             as pl
import oh_module         as ohmd
import copy

from scipy.io.idl        import readsav
from restore             import restore
from mpfit               import mpfit
from scipy.integrate     import quad

## Compute Tex for 1665 line ##
 #
 # params dict data     Data
 # params dict inf408   Info about the Tb_background at 408MHz
 # params bool  xplot   Plot or not
 # params float tau_sig Tau_sigma_limit for Cal. Tex
 # params int  bd       OH665 or OH667
 #
 # return Tex and N(OH)
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def cal_tex_print(data,inf408,xplot=False, tau_sig=2.):
	bg408    = ohmd.read_tbg408_healpy()
 	src_list = list(data.la.srcname)

 	for src in src_list:
	 	n,ell,bee,oh_f1,vlsr1,oh_f2,vlsr2,em_avg1,ab_avg1,em_avg2,ab_avg2 = ohmd.get_src_info(data,src,src_list)
	 	avmin1,avmax1,avmin2,avmax2,evmin1,evmax1,evmin2,evmax2           = ohmd.read_bins_to_cal_bg(n)
	 	xmin, xmax                                                        = ohmd.vel_range(n)

	 	if (src != 'P0428+20'):
			continue
		if (xmin == 0. and xmax == 0.):
			continue

	 	## ABSORPTION and EMISSION DTA: OH1665 & OH1667 ##
	 	xd   = [vlsr1,           vlsr2]
		td   = [ab_avg1,         ab_avg2]
		tde  = [em_avg1,         em_avg2]
		cst  = [3.99757843817,   2.21841824609]
		cst  = [4.26413,         2.36634]
		frq  = [1665.402,        1667.359]
		pfl  = ['../data/gauss_1665_peaks_hiep.txt',   '../data/gauss_1667_peaks.txt']
			
		xd1  = xd[0];  xd2  = xd[1]
		td1  = td[0];  td2  = td[1]
		tde1 = tde[0]; tde2 = tde[1]
		cst1 = cst[0]; cst2 = cst[1]
		frq1 = frq[0]; frq2 = frq[1]
		pfl1 = pfl[0]; pfl2 = pfl[1]
		# xd,td = bin_up(xd,td,nbin=1)

		## CONTIMNUUM BACKGROUND OH65 and OH67 ##
		tbg65 = 2.725 + ohmd.get_tb_408(ell,bee,inf408.tb_408)*(408./frq1)**2.8 # Tbg from 408MHz
		tbg65 = 2.725 + bg408[src]*(408./frq1)**2.8 # Tbg from 408MHz for OH65
		tbg67 = 2.725 + bg408[src]*(408./frq2)**2.8 # Tbg from 408MHz for OH67
		if(src == '3C123'):
			tbg65 = 25.0
			tbg67 = 25.0

		## BASELINES and UNCERTAINTIES ##
		tc65,tc65_er        = ohmd.baseline_from_linear_fit(xd1, td1, avmin1, avmax1, avmin2, avmax2,fit=False)
		bg_off65,bgoff65_er = ohmd.baseline_from_linear_fit(xd1, tde1, evmin1, evmax1, evmin2, evmax2,fit=False)
		trx65               = bg_off65 - tbg65
		tc67,tc67_er        = ohmd.baseline_from_linear_fit(xd2, td2, avmin1, avmax1, avmin2, avmax2,fit=False)
		bg_off67,bgoff67_er = ohmd.baseline_from_linear_fit(xd2, tde2, evmin1, evmax1, evmin2, evmax2,fit=False)
		trx67               = bg_off67 - tbg67

		## 1-SIGMA STANDARD DEVIATION OF Tabspt and Temission ##
		tab65_sigma  = tc65_er
		tem65_sigma  = bgoff65_er
		trx65_sigma  = bgoff65_er
		ton65_sig    = np.sqrt(tem65_sigma**2 + tab65_sigma**2)

		tab67_sigma  = tc67_er
		tem67_sigma  = bgoff67_er
		trx67_sigma  = bgoff67_er
		ton67_sig    = np.sqrt(tem67_sigma**2 + tab67_sigma**2)

		## COMPUTE EXP(-TAU), TAU & 1-SIGMA STANDARD DEVIATION OF TAU ##
		etaud65      = td1/tc65
		taud65       = -np.log(etaud65)
		tau65_sigma  = ohmd.get_tb_sigma(xd1, taud65, avmin1, avmax1, avmin2, avmax2)
		etau65_sigma = np.abs(etaud65)*np.sqrt( (tab65_sigma/td)**2 + (tc65_er/tc65)**2 )
		etau65_sigma = ohmd.get_tb_sigma(xd1, etaud65, avmin1, avmax1, avmin2, avmax2)

		etaud67      = td2/tc67
		taud67       = -np.log(etaud67)
		tau67_sigma  = ohmd.get_tb_sigma(xd2, taud67, avmin1, avmax1, avmin2, avmax2)
		etau67_sigma = np.abs(etaud67)*np.sqrt( (tab67_sigma/td)**2 + (tc67_er/tc67)**2 )
		etau67_sigma = ohmd.get_tb_sigma(xd2, etaud67, avmin1, avmax1, avmin2, avmax2)

		## If PLOT - Print basic Infor of the Source ##
		if(xplot):
			print '***** Basic Infor OH65, OH67*******'
		 	print '	1) Source: '
		 	print '		', src
		 	print '	2) Tcont:'
		 	print '		',tc65,tc67
		 	print '	3) Background of OFF-SOURCE spectrum:'
		 	print '		',bg_off65,bg_off67
		 	print '	4) Radio contimuum obtained from 408MHz:'
		 	print '		',tbg65,tbg67
		 	print '	5) Receiver Temperature, Trx65 = bg_off65 - tbg65, Trx67 = bg_off67 - tbg67:'
		 	print '		',trx65,trx67
		 	print '***** End - Basic Infor *******'
		 	print ''
		
		# VELOCITY-RANGE & INDEXES #
		xmin       = xmin - 4.
		xmax       = xmax + 2.
		xmax65_id  = ohmd.get_vel_index(xd1, xmin)   # xmax_index
		xmin65_id  = ohmd.get_vel_index(xd1, xmax)   # xmin_index
		num_chnl65 = xmax65_id-xmin65_id        # Total number of bins
		vrange65   = [xmin65_id, xmax65_id]
		dv65       = (xmax-xmin)/num_chnl65

		xmax67_id  = ohmd.get_vel_index(xd2, xmin)   # xmax_index
		xmin67_id  = ohmd.get_vel_index(xd2, xmax)   # xmin_index
		num_chnl67 = xmax67_id-xmin67_id        # Total number of bins
		vrange67   = [xmin67_id, xmax67_id]
		dv67       = (xmax-xmin)/num_chnl67

		## (FOR FUN) FIT ABSORPTION LINE FOR TAU, V0 and WIDTH ##
		## Special Case: T0629+10 ##
		guesspar65,base_range65     = ohmd.peak_info(src,pfl1)
		lguess65                    = [tc65] + guesspar65
		x65,etaufit65,etau65,\
		abp65,abper65,npar65,\
		parbase65,pname65,parinfo65 = ohmd.ab_fit(src,xd1,td1,lguess65,xmin65_id,xmax65_id,evmin1,evmax1,evmin2,evmax2)

		guesspar67,base_range67     = ohmd.peak_info(src,pfl2)
		lguess67                    = [tc67] + guesspar67
		x67,etaufit67,etau67,\
		abp67,abper67,npar67,\
		parbase67,pname67,parinfo67 = ohmd.ab_fit(src,xd2,td2,lguess67,xmin67_id,xmax67_id,evmin1,evmax1,evmin2,evmax2)

		## PLOT - Absorption line & Emission line ##
		# if(xplot):
			## Absorption OH65 line ##
			# colors = ['m','g','b','y','c','r','purple','b']
			# plt.plot(x65, etau65 , 'b.-', label='data - Absorption line', ms=10)
			# plt.plot(x65, etaufit65 ,'r-', label='Gaussian fit', lw=2)
			# for i in range(2,len(abp65),3):
			# 	plt.axvline(abp65[i]-abp65[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
			# 	plt.axvline(abp65[i]+abp65[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
			# plt.title(src,fontsize=35)
			# plt.xlabel('$V_{lsr} (km/s)$',fontsize=35)
			# plt.ylabel(r'$e^{-\tau}$',fontsize=35)
			# plt.legend(loc='upper right')
			# plt.grid()
			# plt.show()

			# ## Emission OH65 line ##
			# plt.plot(x65,tde1[xmin65_id:xmax65_id], 'b.-', label='data - Emission line', ms=10)
			# for i in range(2,len(abp65),3):
			# 	plt.axvline(abp65[i]-abp65[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
			# 	plt.axvline(abp65[i]+abp65[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
			# plt.title(src,fontsize=35)
			# plt.xlabel('$V_{lsr} (km/s)$',fontsize=35)
			# plt.ylabel(r'$T (K)$',fontsize=35)
			# plt.legend(loc='upper right')
			# plt.grid()
			# plt.show()

			# ## Absorption OH67 line ##
			# plt.plot(x67, etau67 , 'b.-', label='data - OH67 Absorption line', ms=10)
			# plt.plot(x67, etaufit67 ,'r-', label='Gaussian fit', lw=2)
			# for i in range(2,len(abp67),3):
			# 	plt.axvline(abp67[i]-abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
			# 	plt.axvline(abp67[i]+abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
			# plt.title(src,fontsize=35)
			# plt.xlabel('$V_{lsr} (km/s)$',fontsize=35)
			# plt.ylabel(r'$e^{-\tau}$',fontsize=35)
			# plt.legend(loc='upper right')
			# plt.grid()
			# plt.show()

			# ## Emission OH67 line ##
			# plt.plot(x67,tde2[xmin67_id:xmax67_id], 'b.-', label='data - OH 67 Emission line', ms=10)
			# for i in range(2,len(abp67),3):
			# 	plt.axvline(abp67[i]-abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
			# 	plt.axvline(abp67[i]+abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=3,color='k')
			# plt.title(src,fontsize=35)
			# plt.xlabel('$V_{lsr} (km/s)$',fontsize=35)
			# plt.ylabel(r'$T (K)$',fontsize=35)
			# plt.legend(loc='upper right')
			# plt.grid()
			# plt.show()

		## PARAMS FROM THE FIT: Tau and Width TO CAL. N(OH) ##
		tau65_fit    = []
		v065_fit     = []
		wid65_fit    = []
		taufit65_sig = []
		v0fit65_sig  = []
		widfit65_sig = []
		for i in range(1,len(abp65),3):
			tau65_fit.append(abp65[i])
			v065_fit.append(abp65[i+1])
			wid65_fit.append(abp65[i+2])
			taufit65_sig.append(abper65[i])
			v0fit65_sig.append(abper65[i+1])
			widfit65_sig.append(abper65[i+2])

		tau67_fit    = []
		v067_fit     = []
		wid67_fit    = []
		taufit67_sig = []
		v0fit67_sig  = []
		widfit67_sig = []
		for i in range(1,len(abp67),3):
			tau67_fit.append(abp67[i])
			v067_fit.append(abp67[i+1])
			wid67_fit.append(abp67[i+2])
			taufit67_sig.append(abper67[i])
			v0fit67_sig.append(abper67[i+1])
			widfit67_sig.append(abper67[i+2])

		## CALCULATE Tex, here CHOOSE etau_data ##
		t65_off = tde1                      # Use OFF-source Spectrum
		t65_off = t65_off[xmin65_id:xmax65_id]
		xde65   = xd1[xmin65_id:xmax65_id]
		e_tau65 = etaud65[xmin65_id:xmax65_id]		
		taud65  = taud65[xmin65_id:xmax65_id]

		t67_off = tde2                      # Use OFF-source Spectrum
		t67_off = t67_off[xmin67_id:xmax67_id]
		xde67   = xd2[xmin67_id:xmax67_id]
		e_tau67 = etaud67[xmin67_id:xmax67_id]		
		taud67  = taud67[xmin67_id:xmax67_id]

		toff65  = []
		xe65    = []
		etaue65 = []
		etauf65 = []
		tau65   = []
		for i in range(len(t65_off)):
			if(taud65[i] > tau_sig*tau65_sigma):
				tau65.append(taud65[i])
				toff65.append(t65_off[i])
				xe65.append(xde65[i])
				etaue65.append(e_tau65[i])
				etauf65.append(etaufit65[i])

		toff67  = []
		xe67    = []
		etaue67 = []
		etauf67 = []
		tau67   = []
		for i in range(len(t67_off)):
			if(taud67[i] > tau_sig*tau67_sigma):
				tau67.append(taud67[i])
				toff67.append(t67_off[i])
				xe67.append(xde67[i])
				etaue67.append(e_tau67[i])
				etauf67.append(etaufit67[i])

		toff65  = np.asarray(toff65,  dtype=np.float64)
		xe65    = np.asarray(xe65,    dtype=np.float64)
		etaue65 = np.asarray(etaue65, dtype=np.float64)
		tau65   = np.asarray(tau65,   dtype=np.float64)
		tex65   = ( toff65-trx65-tbg65*etaue65 )/(1.-etaue65)

		toff67  = np.asarray(toff67,  dtype=np.float64)
		xe67    = np.asarray(xe67,    dtype=np.float64)
		etaue67 = np.asarray(etaue67, dtype=np.float64)
		tau67   = np.asarray(tau67,   dtype=np.float64)
		tex67   = ( toff67-trx67-tbg67*etaue67 )/(1.-etaue67)

		# Find the vel-range of the peak from preapred-file #
		peak65        = ohmd.get_peak_vel_range(abp65)
		npeak65       = len(peak65)/2	
		tex65_peak    = [0.]*npeak65
		tex65sig_peak = [0.]*npeak65
		noh65_peak    = [0.]*npeak65
		stau65sig     = [0.]*npeak65
		noh65_sig     = [0.]*npeak65

		peak67        = ohmd.get_peak_vel_range(abp67)
		npeak67       = len(peak67)/2	
		tex67_peak    = [0.]*npeak67
		tex67sig_peak = [0.]*npeak67
		noh67_peak    = [0.]*npeak67
		stau67sig     = [0.]*npeak67
		noh67_sig     = [0.]*npeak67

		ltex65        = {}
		ltau65        = {}
		lstau65_sig   = {}
		for j in range(npeak65):
			ltex65[j]      = []
			ltau65[j]      = []
			lstau65_sig[j] = []

		ltex67        = {}
		ltau67        = {}
		lstau67_sig   = {}
		for j in range(npeak67):
			ltex67[j]      = []
			ltau67[j]      = []
			lstau67_sig[j] = []

		# CAL. Tex FOR EACH PEAK #
		for i in range(0, len(xe65)):
			for k in range(0,len(peak65),2):
				vmin65 = peak65[0+k]
				vmax65 = peak65[1+k]
				if ( (xe65[i]>=vmin65) and (xe65[i]<=vmax65) ) :
					ltex65[k/2].append(tex65[i])
					ltau65[k/2].append(tau65[i])
					lstau65_sig[k/2].append(tau65_sigma**2)

		for i in range(0, len(xe67)):
			for k in range(0,len(peak67),2):
				vmin67 = peak67[0+k]
				vmax67 = peak67[1+k]
				if ( (xe67[i]>=vmin67) and (xe67[i]<=vmax67) ) :
					ltex67[k/2].append(tex67[i])
					ltau67[k/2].append(tau67[i])
					lstau67_sig[k/2].append(tau67_sigma**2)

		for k in range(npeak65):
			tex65_peak[k]    = np.mean(ltex65[k])
			tex65sig_peak[k] = np.std(ltex65[k])
			noh65_peak[k]    = 4.26413*tex65_peak[k]*tau65_fit[k]*wid65_fit[k]
			noh65_sig[k]     = ohmd.cal_noh65_uncertainty(tau65_fit[k], wid65_fit[k], tex65_peak[k], taufit65_sig[k], widfit65_sig[k], tex65sig_peak[k])

		for k in range(npeak67):
			tex67_peak[k]    = np.mean(ltex67[k])
			tex67sig_peak[k] = np.std(ltex67[k])
			noh67_peak[k]    = 2.36634*tex67_peak[k]*tau67_fit[k]*wid67_fit[k]
			noh67_sig[k]     = ohmd.cal_noh67_uncertainty(tau67_fit[k], wid67_fit[k], tex67_peak[k], taufit67_sig[k], widfit67_sig[k], tex67sig_peak[k])

		## PLOT OFF-Source LINE: t_off  = tde ##
		if(xplot):
			# plt.plot(xe65,toff65,'b.-', label='OFF-source spectrum, $T_{on} = T_{abs} + T_{off}$', ms=10)
			# for i in range(0,len(peak65),2):
			# 	plt.axvline(peak65[i], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
			# 	plt.axvline(peak65[i+1], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
			# plt.title(src,fontsize=35)
			# plt.xlabel('$V_{lsr} (km/s)$',fontsize=35)
			# plt.ylabel(r'$T_{on} (K)$',fontsize=35)
			# plt.legend(loc='upper right')
			# plt.grid()
			# plt.show()

			# ## PLOT Tex spectrum ##
			# plt.plot(xe65,tex65, 'r', label='$T_{ex}$', lw=2)
			# plt.plot(xe65,-7000.*np.log(etaue65), 'b', label=r'$7000*e^{-\tau}$', lw=2)
			# for i in range(0,len(peak65),2):
			# 	plt.axvline(peak65[i], ymin=-10., ymax=1000.,linewidth=3,color='k')
			# 	plt.axvline(peak65[i+1], ymin=-10., ymax=1000.,linewidth=3,color='k')
			# plt.ylim(-50.,400.)
			# plt.title(src, fontsize=35 )
			# plt.xlabel('$V_{lsr} (km/s)$', fontsize=35)
			# plt.ylabel('$T_{EX} (K)$',fontsize=35)
			# plt.legend(loc='upper right')
			# plt.grid()
			# plt.show()

		 #    ## CALCULATE Tau ##
			# plt.plot(xe65,-np.log(etaue65), 'b', label=r'${\tau}$', lw=2)
			# for i in range(0,len(peak65),2):
			# 	plt.axvline(peak65[i], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])
			# 	plt.axvline(peak65[i+1], ymin=-10., ymax=1000.,linewidth=1,color=colors[i/2])

			# plt.title(src,fontsize=35)
			# plt.xlabel('$V_{lsr} (km/s)$',fontsize=35)
			# plt.ylabel(r'${\tau}$',fontsize=35)
			# plt.legend(loc='upper right')
			# plt.grid()
			# plt.show()


			## =========================== ##



			## ABS line fit & residuals
			## Plot
			mpl.rcParams['axes.linewidth'] = 2.
			plt.rc('font', weight='bold')
			plt.rc('text', usetex=True)
			plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

			xmin = 0.
			xmax = 14.

			fig1 = plt.figure(figsize=(11.5,12))
			fig1.set_rasterized(True)

			frame1 = fig1.add_axes((.1,.7,.89,.285))

			major_xticks = np.arange(0., 14., 2.)
			minor_xticks = np.arange(0., 14., 1.)
			major_yticks = np.arange(0.984, 1.002, 0.004)
			minor_yticks = np.arange(0.984, 1.002, 0.001)

			plt.plot(x67, etau67, 'k-', linewidth=1, label='')
			plt.plot(x67, etaufit67, 'k-', linewidth=2, label='')

			for i in range(2,len(abp67),3):
				plt.axvline(abp67[i]-abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color='k', ls='--')
				plt.axvline(abp67[i]+abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color='k', ls='--')
			
			frame1.set_xticks(major_xticks)                                                       
			frame1.set_xticks(minor_xticks, minor=True)                                           
			frame1.set_yticks(major_yticks)                                                       
			frame1.set_yticks(minor_yticks, minor=True)
			plt.tick_params(axis='x', labelsize=16, pad=2, labelbottom='off')
			plt.tick_params(axis='y', labelsize=16, pad=2)
			plt.tick_params(which='both', width=2)
			plt.tick_params(which='major', length=6)
			plt.tick_params(which='minor', length=4)

			yy1 = 0.997
			plt.text(0.23, yy1, r'$\mathrm{P0428{+}20}$', fontsize=16, fontweight='bold')
			plt.text(0.23, yy1-0.001, '$\mathrm{(l,b)=(176.81,-18.56)}$', fontsize=14, fontweight='bold')

			plt.text(0.23, yy1-0.0027, '$\mathrm{OH\ (1667)}$', fontsize=14, fontweight='bold')
			plt.text(0.23, yy1-0.0037, r'$\tau\ \ \ \ \ =[0.0015,\ 0.0076]$', fontsize=14, fontweight='bold')
			plt.text(0.23, yy1-0.0049, '$\mathrm{V_{LSR}}\ \ \ \ =[3.60,\ \ 10.70]$', fontsize=14, fontweight='bold')
			plt.text(0.23, yy1-0.0062, '$\mathrm{FWHM}=[1.06,\ \ 1.10]$', fontsize=14, fontweight='bold')
			plt.text(0.23, yy1-0.0073, '$\mathrm{T_{ex}}\ \ \ \ \ \ =[8.424,\ \ 8.420]$', fontsize=14, fontweight='bold')
			
			plt.grid(False)
			plt.ylabel(r'$e^{-\tau}$', fontsize=26)
			plt.xlim([xmin, xmax])
			plt.ylim([0.985, 1.0015])

			## --------------- ##
			frame2=fig1.add_axes((.1,.4,.89,.285))

			major_xticks = np.arange(0., 14., 2.)
			minor_xticks = np.arange(0., 14., 1.)
			major_yticks = np.arange(-0.04, 0.1, 0.02)
			minor_yticks = np.arange(-0.04, 0.1, 0.01)

			plt.plot(x67,tde1[xmin67_id:xmax67_id]-bg_off67+0.165, 'k-', lw=1.5, label='')


			for i in range(2,len(abp67),3):
				plt.axvline(abp67[i]-abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color='k', ls='--')
				plt.axvline(abp67[i]+abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color='k', ls='--')
			

			frame2.set_xticks(major_xticks)                                                       
			frame2.set_xticks(minor_xticks, minor=True)                                           
			frame2.set_yticks(major_yticks)                                                       
			frame2.set_yticks(minor_yticks, minor=True)
			plt.tick_params(axis='x', labelsize=16, pad=2, labelbottom='off')
			plt.tick_params(axis='y', labelsize=16, pad=2)
			plt.tick_params(which='both', width=2)
			plt.tick_params(which='major', length=6)
			plt.tick_params(which='minor', length=4)

			plt.ylabel(r'$\mathrm{T_{exp}\ [K]}$',fontsize=24)
			plt.tick_params(axis='x', labelsize=18, labelbottom='off')
			plt.tick_params(axis='y', labelsize=15)
			# plt.legend(loc='upper left', fontsize=18)
			plt.xlim([xmin, xmax])
			plt.ylim([-0.04, 0.1])



			frame4=fig1.add_axes((.1,.1,.89,.285))

			major_xticks = np.arange(0., 14., 2.)
			minor_xticks = np.arange(0., 14., 1.)
			major_yticks = np.arange(-30., 40., 10.)
			minor_yticks = np.arange(-30., 40., 5.)


			plt.plot(xe67,tex67, 'k', label='', lw=2)
			# plt.plot(xe67,-7000.*np.log(etaue67), 'b', label='', lw=2)
			# for i in range(2,len(abp67),3):
			# 	plt.axvline(abp67[i]-abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color='k', ls='--')
			# 	plt.axvline(abp67[i]+abp67[i+1]/2.,ymin=-10., ymax=1000.,linewidth=2,color='k', ls='--')

			plt.axvspan(0., abp67[2]-abp67[3]/2., ymin=-30., ymax=50., alpha=0.5, color='k')
			plt.axvspan(abp67[2]+abp67[3]/2., abp67[5]-abp67[6]/2., ymin=-30., ymax=50., alpha=0.5, color='k')
			plt.axvspan(abp67[5]+abp67[6]/2., 14., ymin=-30., ymax=50., alpha=0.5, color='k')

			
			# plt.tick_params(axis='y', labelsize=13, pad=7)
			# plt.tick_params(axis='x', labelsize=18)
			# plt.tick_params(which='both', width=1)
			# plt.tick_params(which='major', length=5, top='off')


			frame4.set_xticks(major_xticks)                                                       
			frame4.set_xticks(minor_xticks, minor=True)                                           
			frame4.set_yticks(major_yticks)                                                       
			frame4.set_yticks(minor_yticks, minor=True)
			# plt.tick_params(axis='x', labelsize=16, pad=2, labelbottom='off')
			plt.tick_params(axis='x', labelsize=16, pad=2)
			plt.tick_params(axis='y', labelsize=16, pad=2)
			plt.tick_params(which='both', width=2)
			plt.tick_params(which='major', length=6)
			plt.tick_params(which='minor', length=4)


			plt.ylabel('$\mathrm{T_{ex}\ [K]}$', fontsize=24)
			plt.xlabel('$\mathrm{VLSR\ [km/s]}$', fontsize=24)
			plt.xlim([xmin, xmax])
			plt.ylim([-30., 40.])
			
			# plt.savefig('tex_spec_'+src+'.png', format='png', dpi=400)
			plt.savefig('tex_spec_P0428+20.eps', bbox_inches='tight', pad_inches=0.03, format='eps', dpi=150)
			plt.show()

#============== MAIN ==============#
data   = readsav('../data/makelines.sav')              # data.la
inf408 = readsav('../data/tb_408.sav')                 # l_cntr, b_cntr, tb_408, Continuum at 408MHz
print '================================================================= Results =================================================================='
print 'n 	src 		tau_fit		v0_fit		wid_fit 	tex_peak 	noh_peak 	tau7_fit		v07_fit		wid7_fit 	tex7_peak 	noh7_peak '
print '============================================================================================================================================='
cal_tex_print(data, inf408, xplot=1, tau_sig=0.) # OH65 & OH67

sys.exit()