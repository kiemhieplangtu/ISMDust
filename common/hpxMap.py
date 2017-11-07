from basicImport import *

#######################################################################
# Block for Placnk whole-map data
# Author: Van Hiep
# Version: 08/2017
##################### - START - #######################################


## Cal E(B-V) from Planck 2014 R1.2 #
 #
 # params array tau_map  Map of Tau353
 # params array tauer_map  Map of Tau353_error
 # params array ebv_map  Map of E(B-V)
 # params dict  info     Infor of the sources
 #
 # return list info
 # 
 # version 08/2017
 # Author Van Hiep ##
def cal_ebv_from_planckR12(tau_map, tauer_map, ebv_map, info):
	# Define constants #
	deg2rad = np.pi/180.

	fct     = 1.49e4
	fct_er  = 0.03e4

	## sources
	src    = info['src']  ## src
	xl     = info['l']
	xb     = info['b']
	n      = len(src)

	nside  = hp.get_nside(ebv_map)
	res    = hp.nside2resol(nside, arcmin=False)

	# OK - Go #
	ebv    = np.zeros(n)
	ebver  = np.zeros(n)
	for i in range(n):
		theta = (90.0-xb[i])*deg2rad
		phi   = xl[i]*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (ebv_map[pix] > -0.000001) : # Some pixels not defined
			xtau   = tau_map[pix]
			xtauer = tauer_map[pix]
			val    = ebv_map[pix]

			d1     = fct_er/fct
			d2     = xtauer/xtau
			dr     = np.sqrt(d1**2 + d2**2)
			err    = val*dr

		ebv[i]   = val
		ebver[i] = err

	return ebv, ebver


## Get tau353 and its error #
 #
 # params array tauMap  Map of tau353
 # params array errMap  Error map of tau353
 # params dict  info    Infor of the sources
 #
 # return 1-D arrays tau353 and tauer
 # 
 # version 12/2016
 # Author Van Hiep ##
def get_tau(tauMap, errMap, info):
	# src and l,b
	src = info['src']
	xl  = info['l']
	xb  = info['b']
	n   = len(src)

	# Define constants
	deg2rad = np.pi/180.
	nside   = hp.get_nside(tauMap)
	res     = hp.nside2resol(nside, arcmin=False)

	# OK - Go #
	tau353 = np.zeros(n)
	tauer  = np.zeros(n)
	for i in range(n):
		# Find the values of Tau353 and Err_tau353
		theta = (90.0-xb[i])*deg2rad
		phi   = xl[i]*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if ( (errMap[pix] >= 6.9e-11) and (errMap[pix] <= 0.00081) and (tauMap[pix] > -1.0e30) ): # Checked Invalid error & Some pixels not defined
			tau353[i] = tauMap[pix] 
			tauer[i]  = errMap[pix]
		else:
			print("Error! Error!")
			sys.exit()

	return tau353, tauer


## Get N(HI) from HI4PI #
 #
 # params array HImap  Map of N(HI) from HI4pi
 # params dict  info    Infor of the sources
 #
 # return 1-D arrays N(HI)
 # 
 # version 12/2016
 # Author Van Hiep ##
def get_HI4pi(xmap, info):
	# src and l,b
	src = info['src']
	xl  = info['l']
	xb  = info['b']
	n   = len(src)

	# Define constants
	deg2rad = np.pi/180.
	nside   = hp.get_nside(xmap)
	res     = hp.nside2resol(nside, arcmin=False)

	# OK - Go #
	nhi = np.zeros(n)
	for i in range(n):
		# Find the values of nhi 
		theta = (90.0-xb[i])*deg2rad
		phi   = xl[i]*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if ( xmap[pix] > -1.0e30 ): # Checked Invalid error & Some pixels not defined
			nhi[i] = xmap[pix] 
		else:
			print("Error! Error!")
			sys.exit()

	return nhi/1e20


## Get Radiance #
 #
 # params 1Darray tauMap     Map of tau353
 # params 1Darray tauErrMap  Map of tau353_error
 # params 1Darray rMap       Map of Radiance
 # params dict    info       Infor of the sources
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def get_radiance(tauMap, tauErrMap, rMap, info):
	# src and l,b
	src = info['src']
	xl  = info['l']
	xb  = info['b']
	n   = len(src)

	# Define constants #
	deg2rad = np.pi/180.

	# Constants to cal. error of Radiance from Tau353 (see PLC2014a for details)
	fct     = 0.0276
	fct_er  = 0.00072

	nside   = hp.get_nside(rMap)
	res     = hp.nside2resol(nside, arcmin=False)

	# OK - Go #
	r       = np.zeros(n)
	rer     = np.zeros(n)
	for i in range(n):
		theta = (90.0-xb[i])*deg2rad
		phi   = xl[i]*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (rMap[pix] > -0.000001) : # Some pixels not defined
			xtau   = tauMap[pix]
			xtauer = tauErrMap[pix]
			val    = rMap[pix]

			d1     = fct_er/fct
			d2     = xtauer/xtau
			dr     = np.sqrt(d1**2 + d2**2)
			err    = val*dr           ## Error of Radiance

		r[i]   = val
		rer[i] = err

	return r, rer


## Bin-up the data #
 #
 # params 1Darray  x  (In ascending order)    
 # params 1Darray  xer   
 # params 1Darray  y      
 # params 1Darray  yer   
 # params Int      n  Number of datapoints in a bin   
 #
 # return lists x, xer, y, yer
 # 
 # version 10/2017
 # Author Van Hiep ##
def bin_up_data(x, n, error=False):
	N     = len(x)
	mod   = N%n
	floor = N//n

	chnl  = floor
	if(mod != 0):
		chnl = chnl + 1

	indx = 0
	lst  = range(0,N,n)
	lmax = max(lst)

	xx   = np.zeros(chnl)
	for i in lst:
		nbin = n
		if( (i == lmax) and (mod != 0)) :
			nbin = mod

		if(error):
			xtemp = np.square(x[i:(i+nbin)])
			xtemp = np.sum(xtemp)
			xtemp = np.sqrt(xtemp)/nbin
		else:
			xtemp = np.mean(x[i:(i+nbin)])

		xx[indx] = xtemp
		indx     = indx + 1

	return xx


## Cal N(H) from tau353 #
 #
 # params array tau_map  Map of tau353
 # params array err_map  Error map of tau353
 # params dict  info     Infor of the sources
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def get_nh_from_tau353(tauMap, errMap, info):
	# Define constants #
	deg2rad     = np.pi/180.
	fukui_cf    = 2.10e26
	fk_fact_err = 0.0 #unknown

	## find Planck Conversion Factor (Dust opacity and Its error) ## 
	a, aer  = [1.39249e6, 4.9102e4] #6.6e-27, 0.66e-26, lowNHI
	b, ber  = [0., 0.]              ## If having intercerpt

	# Cal. tau353
	t353, t353er = get_tau(tauMap, errMap, info)

	# Calculate the NH from Planck factor #
	nh   = t353*a+b
	nher = nh_uncert_from_proxies(t353, t353er, a, aer, b, ber, corrEff=0.)

	return nh, nher, t353, t353er

## Cal N(H) from E(B-V) #
 #
 # params 1D array ebv     E(B-V)
 # params 1D array ebver   Error
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def get_nh_from_ebv(ebv, ebver):
	# a,aErr = [5.8e21, 0.0]        ## From Bohlin 1978   
	a, aer  = [113.912, 4.18695]    ## From S&F2011
	b, ber  = [0., 0.]

	nh      = ebv*a+b
	nher    = nh_uncert_from_proxies(ebv, ebver, a, aer, b, ber, corrEff=0.)

	return nh, nher, ebv, ebver

## Cal NH from Radiance #
 #
 # params array r_map    Map of Radiance
 # params dict  info     Infor of the sources
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def get_nh_from_radiance(tauMap, tauErrMap, rMap, info):
	a, aer = [4.65435e11, 0.170979e11]   ## N(H) = a.e31.R + b, NH/e20 = a.e11.R + b
	b, ber = [0., 0.]

	R, Rer = get_radiance(tauMap, tauErrMap, rMap, info)
	R      = R*1e-4     ## radiance W/m2/sr -> w/cm2/sr, xdata*1e11 => 1e7
	Rer    = Rer*1e-4

	nh     = R*a+b
	nher   = nh_uncert_from_proxies(R, Rer, a, aer, b, ber, corrEff=0.)

	return nh, nher, R, Rer

## Get Tb from Haslam #
 #
 # params str map_file  File of Radiance Map
 # params dict  info   Infor of the sources
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def get_tb_haslam(haslam_map, info):
	## sources
	src     = info['src']  ## src

	# Define the width of area #
	beam    = 5.             # Beam = 36'
	dbeam   = beam/120.0     # Beam = 36' -> dbeam = beam/60/2 in degree
	offset  = dbeam          # degree

	nside   = hp.get_nside(haslam_map)
	res     = hp.nside2resol(nside, arcmin=False)
	dd      = res/const.DEG2RAD/10.0

	# OK - Go #
	tbg    = {}
	for i in range(len(src)):
		l  = info['l'][i]
		b  = info['b'][i]

		# Cal. #
		theta = (90.0-b)*const.DEG2RAD
		phi   = l*const.DEG2RAD
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (haslam_map[pix] > -0.000001) : # Some pixels not defined
			val = haslam_map[pix]

		# tbg.append(val)
		temp        = 2.725 + val*(408./1420.4)**2.7 # Tbg from 408MHz 
		tbg[src[i]] = temp

	return tbg

#######################################################################
# End - Block for Placnk whole-map data
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################