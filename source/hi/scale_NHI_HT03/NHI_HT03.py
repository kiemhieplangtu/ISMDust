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

##================= MAIN ========================##

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


### N(HI) ###
ret = {}
for i in range(len(src)):
	woc    = cnm[i]
	tauer  = err_tau[i]
	tser   = err_tspin[i]
	nhier  = 0.   
	if (woc==9):  # Warm component
		tser   = '----'
		woc    = '-'
		tauer  = '----'
		wNHIer = md.get_NWNM_error(tb[i], err_tb[i], del_v[i], err_del_v[i])
		wNHIer = wNHIer*0.01  ## in Unit of 1e20
		string = '{:.2f}   {:.2f}   {:.3f}   {:.3f}   {:.1f}   {:.1f}   {:.2f}   {:.2f}   {:.2f}   {:.2f}   {:.1f}   {:.2f}   {:.4f}   {}   {:.2f}   {:.2f}   {}'\
		.format(tb[i], err_tb[i], tau[i], err_tau[i], v0[i], v0er[i], del_v[i], err_del_v[i], tspin[i], err_tspin[i], tkmax[i], nhi[i], round(wNHIer,4), cnm[i], frac[i], err_frac[i], src[i])
		print string
	else:
		wNHIer = md.get_NCNM_error(tau[i], tauer, tspin[i], err_tspin[i], del_v[i], err_del_v[i])
		wNHIer = wNHIer*0.01  ## in Unit of 1e20

		string = '{:.2f}   {:.2f}   {:.3f}   {:.3f}   {:.1f}   {:.1f}   {:.2f}   {:.2f}   {:.2f}   {:.2f}   {:.1f}   {:.2f}   {:.4f}   {}   {:.2f}   {:.2f}   {}'\
		.format( tb[i], err_tb[i], tau[i], err_tau[i], v0[i], v0er[i], del_v[i], err_del_v[i], tspin[i], err_tspin[i], tkmax[i], nhi[i], round(wNHIer,4), cnm[i], frac[i], err_frac[i], src[i] )
		print string

	if src[i] not in ret.keys():
		ret[src[i]] = {}

		ret[src[i]]['nhi']   = nhi[i]
		ret[src[i]]['nhier'] = wNHIer**2
	else:
		ret[src[i]]['nhi']   = ret[src[i]]['nhi']   + nhi[i]
		ret[src[i]]['nhier'] = ret[src[i]]['nhier'] + wNHIer**2


### Cal. Error - Squared-root ##
for sc in ret:
	ret[sc]['nhier'] = np.sqrt(ret[sc]['nhier'])

## print to: fit_params.dat ##