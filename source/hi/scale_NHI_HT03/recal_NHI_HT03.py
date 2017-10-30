import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md

from scipy.io.idl        import readsav
from restore             import restore


## Get infor of 79 los #
 #
 # params str fname  File-name
 # params dict  info   Infor of the sources
 #
 # return list info
 # 
 # version 09/2017
 # Author Van Hiep ##
def read_79_info(fname='../data/79src_radec_lb.txt'):
	cols  = ['l','b','src','ra','dec']
	fmt   = ['f','f','s',  'f',  'f' ]
	dat   = restore(fname, 3, cols, fmt)
	inf   = dat.read(asarray=True)

	return inf


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
	# Define constants #
	deg2rad = np.pi/180.

	## sources
	src     = info['src']  ## src

	# Define the width of area #
	beam    = 5.             # Beam = 36'
	dbeam   = beam/120.0     # Beam = 36' -> dbeam = beam/60/2 in degree
	offset  = dbeam          # degree

	nside   = hp.get_nside(haslam_map)
	res     = hp.nside2resol(nside, arcmin=False)
	dd      = res/deg2rad/10.0

	# OK - Go #
	tbg    = {}
	for i in range(len(src)):
		l  = info['l'][i]
		b  = info['b'][i]

		# Cal. #
		theta = (90.0-b)*deg2rad
		phi   = l*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (haslam_map[pix] > -0.000001) : # Some pixels not defined
			val = haslam_map[pix]

		# tbg.append(val)
		temp        = 2.725 + val*(408./1420.4)**2.7 # Tbg from 408MHz 
		tbg[src[i]] = temp

	return tbg


##================= MAIN ========================##
map_file = os.getenv("HOME")+'/hdata/oh/haslam408_dsds_Remazeilles2014.fits'
factor   = (408./1420.4)**2.8
tb408    = hp.read_map(map_file,field=0, nest=False, hdu=1, h=False, verbose=False)
info79   = read_79_info(fname='../data/79src_radec_lb.txt')
tbg21    = get_tb_haslam(tb408, info79)

inf408   = readsav('../../oh/data/tb_408.sav')


## Read fit params of HI components from Carl papers ##
cols  = ['t_peak','err_t_peak','tau','err_tau','v0','err_v0','del_v','err_del_v', 'tspin', 'err_tspin', 'tkmax', 'nhi', 'cnm', 'frac', 'err_frac', 'src']
fmt   = ['f',     'f',          'f',  'f',      'f', 'f',     'f',    'f',         'f',      'f',         'f',     'f',   'i',    'f',   'f',        's']
dat   = restore('../result/component_fit_params.txt', 3, cols, fmt)
inf   = dat.read(asarray=True)

src       = inf['src']
frac      = inf['frac']
cnm       = inf['cnm']
nhi       = inf['nhi']
err_tau   = inf['err_tau']
err_tspin = inf['err_tspin']
v0        = inf['v0']
v0er      = inf['err_v0']
tkmax     = inf['tkmax']
err_frac  = inf['err_frac']

del_v     = inf['del_v']
err_del_v = inf['err_del_v']
tau       = inf['tau']
tspin     = inf['tspin']
tb        = inf['t_peak']
err_tb    = inf['err_t_peak']


### Re-calc. N(HI) ###
A   = 1.26
ret = {}
for i in range(len(src)):
	newNHI = 0.
	woc    = cnm[i]
	tauer  = err_tau[i]
	tser   = err_tspin[i]
	nhier  = 0.   
	if (woc==9):  # Warm component
		tser   = '----'
		woc    = '-'
		tauer  = '----'
		newNHI = A*nhi[i]  ## already in 1e20
		wNHIer = A*md.get_NWNM_error(tb[i], err_tb[i], del_v[i], err_del_v[i])
		wNHIer = wNHIer*0.01  ## in Unit of 1e20
		string = '{:.2f}   {:.2f}   {:.3f}   {:.3f}   {:.1f}   {:.1f}   {:.2f}   {:.2f}   {:.4f}   {:.4f}   {:.1f}   {:.4f}   {:.4f}   {}   {:.2f}   {:.2f}   {}'\
		.format(tb[i], err_tb[i], tau[i], err_tau[i], v0[i], v0er[i], del_v[i], err_del_v[i], tspin[i], err_tspin[i], tkmax[i], newNHI, round(wNHIer,4), cnm[i], frac[i], err_frac[i], src[i])
		print string
	else:
		etau   = np.exp(-tau[i])
		tbi    = tbg21[src[i]]

		newTs  = A*tspin[i] + (A-1.0)*tbi*etau/(1.0-etau)
		wTser  = md.new_Tspin_error(tbi, A, tser, tau[i], tauer)
		newNHI = 1.93988*tau[i]*newTs*del_v[i]
		newNHI = 0.01*newNHI   ## in Unit of 1e20
		wNHIer = md.get_NCNM_error(tau[i], tauer, newTs, wTser, del_v[i], err_del_v[i])
		wNHIer = wNHIer*0.01  ## in Unit of 1e20

		string = '{:.2f}   {:.2f}   {:.3f}   {:.3f}   {:.1f}   {:.1f}   {:.2f}   {:.2f}   {:.4f}   {:.4f}   {:.1f}   {:.4f}   {:.4f}   {}   {:.2f}   {:.2f}   {}'\
		.format( tb[i], err_tb[i], tau[i], err_tau[i], v0[i], v0er[i], del_v[i], err_del_v[i], round(newTs,4), round(wTser,4), tkmax[i], round(newNHI,4), round(wNHIer,4), cnm[i], frac[i], err_frac[i], src[i] )
		print string

	if src[i] not in ret.keys():
		ret[src[i]] = {}

		ret[src[i]]['nhi']   = newNHI
		ret[src[i]]['nhier'] = wNHIer**2
	else:
		ret[src[i]]['nhi']   = ret[src[i]]['nhi']   + newNHI
		ret[src[i]]['nhier'] = ret[src[i]]['nhier'] + wNHIer**2


### Cal. Error - Squared-root ##
for sc in ret:
	ret[sc]['nhier'] = np.sqrt(ret[sc]['nhier'])

## print to: fit_params_recal.dat ##