import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import module            as md

from restore             import restore

## Cal NH from Radiance #
 #
 # params array tauMap    Map of Tau353
 # params array tauerMap  Map of Tau353_error
 # params array rMap      Map of Radiance
 # params dict  info      Information of 94 sources
 # params dict  avInf     Information of E(B-V) and Av for 21 OH sources
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def nh_from_radiance(tauMap, tauerMap, rMap, info, avInf):
	## sources
	src   = info['src']
	nhi   = info['nhi']
	nhier = info['nhi_er']
	cnm   = info['cnm']
	cnmer = info['cnm_er']

	# Define constants #
	deg2rad  = np.pi/180.
	fct      = 0.0276   ## Tau and Radiance
	fct_er   = 0.00072

	## Radiance to N(H)
	cf,cf_er = [4.05e31, 0.32e31]   ## From Radiance vs N(HI), 19 src without CO & OH and 16 src with low N(HI)
	of,of_er = [0.09e20, 0.21e20]

	nside    = hp.get_nside(rMap)
	res      = hp.nside2resol(512, arcmin=True)
	dd       = res/deg2rad/10.0

	print "Map resolution: ", res

	# OK - Go #
	rd       = []
	rder     = []

	rnh2     = []
	rnh2_er  = []

	rnh      = []
	rnh_er   = []

	rnhi     = []
	rnhi_er  = []

	rav      = []
	rav_er   = []

	rsrc     = []

	rcnm     = []
	rcnm_er  = []

	xl       = []
	xb       = []
	for i in range(len(src)):
		# Find the values of ebv353 and Err_ebv353 in small area #
		l     = info['l'][i]
		b     = info['b'][i]

		# Cal. #
		theta = (90.0-b)*deg2rad
		phi   = l*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (rMap[pix] > -0.000001) : # Some pixels not defined
			xtau   = tauMap[pix]
			xtauer = tauerMap[pix]
			rad_i  = rMap[pix]

			d1     = fct_er/fct
			d2     = xtauer/xtau
			dr     = np.sqrt(d1**2 + d2**2)
			err    = rad_i*dr
			err    = xtau*fct*dr

		rad_i = 1e-4 * rad_i  ## radiance W/m2/sr -> w/cm2/sr, xdata*1e11 => 1e7
		err   = 1e-4 * err  ## radiance W/m2/sr -> w/cm2/sr, xdata*1e11 => 1e7

		rd.append(rad_i*1e11) ## [1e11 w/cm2/sr]
		rder.append(err*1e11)
	   
		## Calculate the NH from Radiance #
		n_h    = cf*rad_i + of
		n_h    = n_h/1e20
		nh_er  = md.uncertainty_of_product(cf, rad_i, cf_er, err)
		nh_er  = md.nh_uncert_from_proxies(cf, rad_i, cf_er, err)
		nh_er  = np.sqrt(nh_er**2 + of_er**2)
		nh_er  = nh_er/1e20

		## N(H2) = (NH-NHI)/2 ##
		nh2    = (n_h-nhi[i])/2.
		nh2_er = 0.5*md.uncertainty_of_diff(nh_er, nhier[i])

		string = '{:10s} {:08.4f}   {:08.4f}   {:10.6f}   {:10.6f}   {:08.4f}   {:08.4f}   {:08.4f}   {:08.4f}   {:08.4f}   {:08.4f}'\
		.format(src[i], l, b, rad_i*1e11, err*1e11, n_h, nh_er, nhi[i], nhier[i], nh2, nh2_er)
		print string

		## N(H2) ##
		rnh2.append(nh2)
		rnh2_er.append(nh2_er)

		## N(H) ##
		rnh.append(n_h)
		rnh_er.append(nh_er)

		## N(HI) ##
		rnhi.append(nhi[i])
		rnhi_er.append(nhier[i])

		## CNM ##
		rcnm.append(cnm[i])
		rcnm_er.append(cnmer[i])

		## l,b ##
		xl.append(l)
		xb.append(b)

	# md.write2file('nh_from_rad_94src.txt', Str)
	return src, xl, xb, rd, rder, rnh, rnh_er, rnhi, rnhi_er, rnh2, rnh2_er, rcnm, rcnm_er

## Get N(H) from Radiance #
 #
 # params str mapFile File of maps
 # params dict info   Information of sources
 # params dict avInf  Information of E(B-V) and Av for 21 OH sources
 #
 # return void
 #
 # version 08/2017 
 # Author Van Hiep
 ##	
def nh_from_rad(mapFile, info94, avInf):
	src   = info94['src']
	nhi   = info94['nhi']
	nhier = info94['nhi_er']
	xl    = info94['l']
	xb    = info94['b']

	## Radiance map R1.2 ##
	tauMap                        = hp.read_map(mapFile, field = 0)
	tauerMap                      = hp.read_map(mapFile, field = 1)
	rMap                          = hp.read_map(mapFile, field = 3)
	print 'src l   b   radiance[1e11]   errRad[1e11]   N(H)   err   N(HI)   err   N(H2)   err'
	src, xl, xb, rd, rder, \
	rnh, rnh_er, rnhi, rnhi_er, \
	rnh2, rnh2_er, rcnm, rcnm_er  = nh_from_radiance(tauMap, tauerMap, rMap, info94, avInf)

	print len(rnh2)
	print len(rnh2_er)

	## NH vs Radiance ##	
	plt.errorbar(rd, rnh, xerr=rder, yerr=rnh_er, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
	plt.title('N$_{H}$ (from Hiep) vs R', fontsize=30)
	plt.xlabel('Rad', fontsize=35)
	plt.ylabel('$N_{H}$', fontsize=35)
	plt.grid(True)
	# plt.ylim(-5e-7, 1e-6)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=15)
	plt.show()


##================= MAIN ========================##
## Filename of the map
pth     = os.getenv("HOME")+'/hdata/dust/'
mapFile = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'  ## E(B-V) from Planck r.12, IRAS ~5'

avInf   = md.read_av_for_oh_src(fname = '../ebv2nh/data/ebv_sfd98_sf2011_for_oh_src.txt')
info94  = md.read_nhi_93src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_scaled.txt')

## cal N(H)
nh_from_rad(mapFile, info94, avInf)