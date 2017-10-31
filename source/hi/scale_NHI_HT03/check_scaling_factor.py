import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

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

	# xl    = inf['l']
	# xb    = inf['b']
	# src   = inf['src']
	# ra    = inf['ra']
	# dec   = inf['dec']

	return inf

## Get infor of 16 common los #
 #
 # params str fname  File-name
 # params list  info   Infor of the sources
 #
 # return list info
 # 
 # version 09/2017
 # Author Van Hiep ##
def read_16common_src(fname='../sponge/16common_src.txt'):
	cols  = ['src', 'l','b', 'msnhi', 'msnhier', 'spnhi', 'spnhier']
	fmt   = ['s',   'f', 'f',  'f',    'f',       'f',      'f'    ]
	dat   = restore(fname, 2, cols, fmt)
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

## Read 16 common src ##
cdat     = read_16common_src(fname='../sponge/old/16common_src.txt')
csrc     = cdat['src']
csrc     = csrc.tolist()
msnhi    = cdat['msnhi']
msnhier  = cdat['msnhier']
spnhi    = cdat['spnhi']
spnhier  = cdat['spnhier']
xl       = cdat['l']
xb       = cdat['b']


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

del_v     = inf['del_v']
err_del_v = inf['err_del_v']
tau       = inf['tau']
tspin     = inf['tspin']
tb        = inf['t_peak']
err_tb    = inf['err_t_peak']


# 	# print src[i], woc, tau[i], tspin[i], xTs
# 	string = '{}   {:10s}   {}   {}   {:.2f}   {:.4f}   {:.4f}'\
# 	.format(i,    src[i],   woc, tau[i], tspin[i], xTs, xTs/tspin[i]   )
# 	print string


### Re-calc. N(HI)
A   = 1.26
ret = {}
for i in range(len(src)):
	if (src[i] in csrc):
		sidx   = csrc.index(src[i]) 
		newNHI = 0.
		woc    = cnm[i]
		tauer  = err_tau[i]
		tser   = err_tspin[i]
		nhier  = 0.   
		if (woc==9):
			tser   = '----'
			woc    = '-'
			tauer  = '----'
			newNHI = A*nhi[i]
			wNHIer = A*md.get_NWNM_error(tb[i], err_tb[i], del_v[i], err_del_v[i])
			wNHIer = wNHIer*0.01  ## in Unit of 1e20
		else:
			etau   = np.exp(-tau[i])
			tbi    = tbg21[src[i]]
			# tbi    = 2.725+md.get_tb_408(xl[sidx],xb[sidx],inf408.tb_408)*(408./1420.4)**2.7 # Tbg from 408MHz

			newTs  = A*tspin[i] + (A-1.0)*tbi*etau/(1.0-etau)
			wTser  = md.new_Tspin_error(tbi, A, tser, tau[i], tauer)
			newNHI = 0.01*1.93988*tau[i]*newTs*del_v[i]
			wNHIer = 0.01*md.get_NCNM_error(tau[i], tauer, newTs, wTser, del_v[i], err_del_v[i])
			print tspin[i], newTs, newTs/tspin[i]

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

newcl   = np.zeros(16)
newcler = np.zeros(16)
for i in range(len(csrc)):
	sc         = csrc[i]
	newcl[i]   = ret[sc]['nhi']
	newcler[i] = ret[sc]['nhier']




### Fit ###

xfit, yfit, mu, sig, m, ea, b, eb         = md.do_linODRfit(msnhi, newcl, msnhier, newcler, lguess=[1.25, 0.5], plot=False)
xfit0, yfit0, mu0, sig0, m0, ea0, b0, eb0 = md.do_linODRfit(msnhi, spnhi, msnhier, spnhier, lguess=[1.25, 0.5])


### Plot ###
plt.figure(1, figsize=(16, 14))
plotboth = 1
if(plotboth):
	plt.plot(xfit, mu, '-r', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

	plt.plot(xfit0, mu0, '-k', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
	plt.fill_between(xfit0, mu0 - sig0, mu0 + sig0, color='0.5', alpha=0.5)

	plt.errorbar(msnhi, newcl, xerr=msnhier, yerr=newcler, color='r', marker='o', ls='', markersize=6, markeredgecolor='r', markeredgewidth=1, label='data')
	plt.errorbar(msnhi, spnhi, xerr=msnhier, yerr=spnhier, color='k', marker='o', ls='', markersize=6, markeredgecolor='k', markeredgewidth=1, label='data')

	plt.plot([0,60],[0,60], 'k--', label='$x=y$')
	plt.show()
else:
	plt.plot(xfit, mu, '-r', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
	plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

	plt.errorbar(msnhi, newcl, xerr=msnhier, yerr=newcler, color='r', marker='o', ls='', markersize=6, markeredgecolor='r', markeredgewidth=1, label='data')

	plt.plot([0,60],[0,60], 'k--', label='$x=y$')
	plt.show()