import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

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

## Read info of XOH from tau, Ebv, Radiance #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 4/2017
 # Author Van Hiep ##
def read_CO_spec(fname = '3c18f_bsl_12co10.dat'):
	cols = ['v','T']
	fmt  = ['f', 'f']
	data = restore(fname, 0, cols, fmt)
	dat  = data.read(asarray=True)
	v    = dat['v']
	T    = dat['T']
	return [v,T]


##================= MAIN ========================##
src = ['3c18', '3c33', '3c142.1', '3c138', '3c409', '3c78', '3c132', '3c310', '3c315', '3c109', 
        '3c237', '3c318', '3c64', '3c454.0', 'p0531+19', '3c333', '3c192', '3c98', '3c273', 'p1055+20', '3c433',
        '3c120', '3c245', '3c348', '4c13.65', '3c274.1', '3c207']
## data Files
pth = os.getenv("HOME")+'/program/class/newly_reduced/'

for sc in src:
	n1  = sc+'f_bsl_12co10.dat'
	n2  = sc+'f_bsl_13co10.dat'
	n3  = sc+'f_bsl_c18o10.dat'
	f1  = pth + n1
	f2  = pth + n2
	f3  = pth + n3

	v1, T1 = read_CO_spec(f1)
	v2, T2 = read_CO_spec(f2)
	v3, T3 = read_CO_spec(f3)

	vlim1  = -12.
	vlim2  = -5.
	rms1, v_fltr1, T_fltr1 = md.get_1sigma_spec(v1,T1, vlim1, vlim2, vmin=-50., vmax=50.)
	rms2, v_fltr2, T_fltr2 = md.get_1sigma_spec(v2,T2, vlim1, vlim2, vmin=-50., vmax=50.)
	rms3, v_fltr3, T_fltr3 = md.get_1sigma_spec(v3,T3, vlim1, vlim2, vmin=-50., vmax=50.)

	plt.plot(v1, T1, 'r-', lw=1, label='')
	# plt.plot(v_fltr1, T_fltr1,'b-')
	plt.axhline(y=2.*rms1, lw=2)
	plt.text(-49., -0.2, '12CO(1-0)', color='b',fontsize=14)

	plt.plot(v2, T2+1., 'r-', lw=1, label='')
	# plt.plot(v_fltr2, T_fltr2+1.,'b-')
	plt.axhline(y=2.*rms2+1., lw=2)
	plt.text(-49., 0.8, '13CO(1-0)', color='b',fontsize=14)


	plt.plot(v3, T3+2., 'r-', lw=1, label='')
	# plt.plot(v_fltr3, T_fltr3+2.,'b-')
	plt.axhline(y=2.*rms2+2., lw=2)
	plt.text(-49., 1.8, 'C18O(1-0)', color='b',fontsize=14)

	plt.grid()
	plt.title(sc, fontsize = 35)
	plt.ylabel('$T_{b} (K)$', fontsize = 35)
	plt.xlabel('VLSR (km/s)', fontsize = 35)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)
	plt.xlim(-50., 50.)
	plt.ylim(-1., 3.)

	# plt.legend(loc='upper right', fontsize = 18)
	# plt.axvline(x=60., lw=4)
	plt.savefig(sc+'.png', format='png', dpi=600)
	plt.show()


