import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md
import copy

from numpy               import array
from restore             import restore
from plotting            import cplot
from mpfit               import mpfit
from scipy.odr           import *

## Linear fucntion ##
 #
 # params list/float p Parameters
 # params 
 #
 # return 
 #
 # version 03/2017 
 # author Nguyen Van Hiep ##
def lin_fc(p, x):
     m = p
     return m*x

## Linear ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
# def tb_exp(x,y,bg,tau,v0,wid,tex,cont,err=1.):
def myfunc(p, fjac=None, x=None, y=None, err=None):
	# model  = p[0] * x + p[1]
	model  = p[0] * x
	status = 0
	return [status, (y - model) / err]

## Read info of OH sources #
 # l,b, noh, noh_er
 #
 # params string fname Filename
 # return dict info
 # 
 # version 4/2017
 # Author Van Hiep ##
def read_noh(fname = '../../oh/result/total_noh65_21src.txt'):
	cols = ['idx','src','l', 'b', 'noh', 'noh_er']
	fmt  = ['i', 's',   'f', 'f', 'f',   'f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	noh  = dat['noh']
	er   = dat['noh_er']
	src  = dat['src']
	l    = dat['l']
	b    = dat['b']

	ret  = {}
	for i in range(0,len(src)):
		# if dat['src'][i] not in ret.keys():
		ret[src[i]] = {}
		ret[src[i]]['noh']   = noh[i]
		ret[src[i]]['l']     = l[i]
		ret[src[i]]['b']     = b[i]
		ret[src[i]]['noher'] = er[i]

	return ret

## Cal NH from Radiance #
 #
 # params array r_map    Map of Radiance
 # params dict  info     Infor of the sources
 # params dict  noh      Infor of OH sources
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def nh_from_radiance(tau_map, tauer_map, r_map, info, noh):
	## sources
	src   = info['src']
	nhi   = info['nhi']
	nhier = info['nhi_er']

	ohsc  = []
	for sc in noh:
		if(sc in src):
			ohsc.append(sc)

	# Define constants #
	deg2rad  = np.pi/180.
	fct      = 0.0276   ## Tau and Radiance
	fct_er   = 0.00072

	# Define the width of area #
	beam     = 5.            # Beam = 36'
	dbeam    = beam/120.0     # Beam = 36' -> dbeam = beam/60/2 in degree
	offset   = dbeam          # degree

	nside    = hp.get_nside(r_map)
	res      = hp.nside2resol(nside, arcmin=False)
	dd       = res/deg2rad/10.0

	cf,cf_er = [4.05e31, 0.32e31]   ## From Radiance vs N(HI), 26 src without CO and 21 src with low N(HI)
	of,of_er = [0.09e20, 0.21e20]

	# OK - Go #
	rd       = []
	rder     = []

	nh       = []
	nher     = []

	roh      = []
	roh_er   = []
	rnh2     = []
	rnh2_er  = []
	rsrc     = []
	for i in range(0, len(src)):
		if (src[i] not in ohsc):
			continue

		# Find the values of ebv353 and Err_ebv353 in small area #
		ci_i  = []

		l     = info['l'][i]
		b     = info['b'][i]

		# Cal. #
		theta = (90.0-b)*deg2rad
		phi   = l*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (r_map[pix] > -0.000001) : # Some pixels not defined
			xtau   = tau_map[pix]
			xtauer = tauer_map[pix]
			val    = r_map[pix]

			d1     = fct_er/fct
			d2     = xtauer/xtau
			dr     = np.sqrt(d1**2 + d2**2)
			err    = val*dr

		val = 1e-4 * val  ## radiance W/m2/sr -> w/cm2/sr, xdata*1e11 => 1e7
		err = 1e-4 * err  ## radiance W/m2/sr -> w/cm2/sr, xdata*1e11 => 1e7

		rd.append(val)
		rder.append(err)
	   
		## Calculate the NH from E(B-V) #
		n_h   = cf*val + of
		nh_er = md.uncertainty_of_product(cf, val, cf_er, err)
		nh_er = np.sqrt(nh_er**2 + of_er**2)

		print src[i],val, n_h, nhi[i], nh_er, nhier[i]

		## N(H2) = (NH-NHI)/2 ##
		nh2   = (n_h-nhi[i]*1e20)/2.

		## 3sources with NHI > NH
		if(nh2 < 0.):
			continue

		nh2_er = 0.5*md.uncertainty_of_diff(nh_er, nhier[i]*1e20)

		## N(H2) ##
		rnh2.append(nh2)
		rnh2_er.append(nh2_er)

		## N(OH) ##
		roh.append( noh[src[i]]['noh'] )
		roh_er.append( noh[src[i]]['noher'] )
		rsrc.append( src[i] )

		# print src[i], n_h, nh_er, nhi[i], nhier[i], nh2, nh2_er, noh[src[i]]['noh'], noh[src[i]]['noher']

	return rnh2, rnh2_er, roh, roh_er, rsrc

## Get tau353 values and err_tau353 values #
 #
 # params str map_file File of maps
 # params dict info   Information of sources
 # params dict lownhi Information of 21 OH sources
 #
 # return void
 #
 # version 11/2016 
 # Author Van Hiep
 ##	
def noh_vs_nh2(map_file, info, noh):
	src   = info['src']
	nhi   = info['nhi']
	# oh    = info['oh']
	nhier = info['nhi_er']
	xl    = info['l']
	xb    = info['b']

	## Radiance map R1.2 ##
	tau_map                        = hp.read_map(map_file, field = 0)
	tauer_map                      = hp.read_map(map_file, field = 1)
	r_map                          = hp.read_map(map_file, field = 3)
	nh2, nh2_er, roh, roh_er, rsrc = nh_from_radiance(tau_map, tauer_map, r_map, info, noh)

	# sys.exit()

	print len(nh2)
	print len(nh2_er)
	print len(roh)
	print len(roh_er)

	## To Plot ##
	xdata = roh
	ydata = nh2

	# Error bar for x-axis and y-axis
	xerr  = roh_er
	yerr  = nh2_er

	########### MPFIT ############
	xdata = 1.0e14*np.array(xdata)
	ydata = np.array(ydata)

	# Error bar for x-axis and y-axis
	xerr  = 1.0e14*np.array(xerr)
	yerr  = np.array(yerr)

	## Histogram ##
	xoh   = xdata/ydata
	print 'Mean: ', xoh.mean()
	for x in xoh:
		print 3, x*1e8 ## xoh_from_radiance.txt
	
	for i in range(len(nh2)):
		print i, rsrc[i], roh[i], nh2[i], xoh[i]

	xmin  = xoh.min()
	xmax  = xoh.max()
	plt.hist(xoh, alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(xmin,xmax), bins=60, lw=3)

	plt.title('', fontsize=22)
	plt.ylabel('Number', fontsize=22,fontweight='bold')
	plt.xlabel('$X_{OH} = N(OH)/N(H_{2})$', fontsize=22, fontweight='bold')

	plt.tick_params(axis='x', labelsize=22, pad=8)
	plt.tick_params(axis='y', labelsize=22)
	plt.tick_params(which='both', width=2)
	plt.tick_params(which='major', length=9)
	plt.tick_params(which='minor', length=4)
	plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))


	plt.legend(loc='upper right', fontsize=22)
	plt.grid(False)

	plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(HI) (<3.0e20)', color='blue', fontsize=20)

	plt.tight_layout()
	plt.show()



##================= MAIN ========================##
## Filename of the map
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'  ## E(B-V) from Planck r.12, IRAS ~5'

# Info of 26 sources with no CO - l/b/name && 23 src low NHI #
# info     = read_nhi_94src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')
info     = md.read_nhi_93src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_scaled.txt')

noh      = read_noh(fname = '../../oh/result/total_noh65_21src.txt')
# noh    = read_noh(fname = '../../oh/result/total_noh67_21src.txt')

print noh

## cal N(H2)
noh_vs_nh2(map_file, info, noh)