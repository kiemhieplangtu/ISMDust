import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class
import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl

import operator
from restore             import restore
from mpl_toolkits.axes_grid1 import make_axes_locatable

## Read NHI from 94src #
 #
 # params string fname Filename
 # return dict info of N(HI)
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_nhi_94src(fname = '../../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt'):
	cols = ['indx', 'src', 'l', 'b', 'nhi', 'nhi_er', 'thin', 'thin_er', 'cnm', 'cnm_er', 'wnm', 'wnm_er']
	fmt  = ['i',     's',   'f','f',  'f',   'f',     'f',     'f',        'f',   'f',      'f',   'f'   ]
	data = restore(fname, 2, cols, fmt)
	return data.read()


## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict infocd 
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_23src_with_OH(fname = '../../../oh/result/23src_with_oh.txt'):
	cols = ['src','l','b']
	fmt  = ['s',  'f','f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()

	src  = dat['src']
	xl   = dat['l']
	xb   = dat['b']

	print xl

	for i in range(len(xl)):
		if (xl[i]>180):
			xl[i] = xl[i]-360.

	print xl
	return src, xl, xb
##================= MAIN ========================##
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
sfd_map  = pth + 'lambda_sfd_ebv.fits'
quantity = 'tau353'
map_unit = ''

# map_file = 'data/HFI_SkyMap_353_2048_R2.02_full.fits'
# quantity = 'w353'
# map_unit = 'K_CMB'


# info  = read_info_no_co('26src_no_co_info.dat')
info = read_nhi_94src(fname = '../../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')
src  = info['src']

# Define the width of area #
beam   = 5.0            # Beam = 5' = map_resolution
dbeam  = beam/120.0     # Beam = 3.5' -> dbeam = beam/60/2 in degree
offset = dbeam          # degree

# Define constants #
nside     = 2048
deg2rad   = np.pi/180.

tau_map   = hp.read_map(map_file, field = 0)
# ebv_map   = hp.read_map(map_file, field = 2)
ebv_map   = hp.read_map(sfd_map, field = 0)
rad_map   = hp.read_map(map_file, field = 3)
tmp_map   = hp.read_map(map_file, field = 4)
nside     = hp.get_nside(tau_map)
res       = hp.nside2resol(nside, arcmin=False)
dd        = res/deg2rad/15.0

# ohsrc, xl, xb = read_23src_with_OH(fname = '../../../oh/result/23src_with_oh_5.dat')
ohsrc, xl, xb = read_23src_with_OH(fname = '../../../oh/result/4src_noCO_noOH_for_plotting.txt')

dg = 1.5
k  = 0
sc = ['3C18', '3C105', '3C109', '3C75']
tt = [r'$\tau_{353}$' + '\n'+ '$[10^{-6}]$', 'E(B-V)'+ '\n'+ '$[10^{0}]mag$', 'Radiance'+ '\n'+ '$[10^{-8}]Wm^{-2}sr^{-1}$', 'Temperature'+ '\n'+ '$K$']
tt = [r'$\tau_{353}$', 'E(B-V)', 'Radiance', 'Temperature']

## Plot
fig = plt.figure(figsize=(13,10))
fig.set_rasterized(True)
for i in range(len(ohsrc)):
	if k > 0:
		tt =['','','','']

	if(i>3):
		continue

	cbr = True

	l  = xl[i]
	b  = xb[i]

	lonra1 = xl[i]-dg
	lonra2 = xl[i]+dg
	if(lonra1<-180.):
		lonra1 = -180.
	if(lonra2>180.):
		lonra2=180.

	lonra = [ lonra1, lonra2 ]
	latra = [ xb[i]-dg, xb[i]+dg ]

	hp.cartview(tau_map*1e6, title=tt[0], coord='G', unit='$(10^{-6})$', norm=None, xsize=800, lonra=lonra, latra=latra, sub=(4,4,k+1),\
	notext = True, cbar=cbr, margins=(0.005,0.00001,0.005,0.005) )
	# hp.projtext(ll, b, '  (' + str(round(ll,2)) + ', ' + str(round(b,2)) + ')', lonlat=True, coord='G', fontsize=15, weight='bold', color='y')
	hp.cartview(ebv_map, title=tt[1], coord='G', unit='$mag$', norm=None, xsize=800, lonra=lonra, latra=latra, sub=(4,4,k+2),\
	notext = True, cbar=cbr, margins=(0.005,0.00001,0.005,0.005))
	hp.cartview(rad_map*1e8, title=tt[2], coord='G', unit='$(10^{-8})Wm^{-2}sr^{-1}$', norm=None, xsize=800, lonra=lonra, latra=latra, sub=(4,4,k+3),\
	notext = True, cbar=cbr, margins=(0.005,0.00001,0.005,0.005) )
	hp.cartview(tmp_map, title=tt[3], coord='G', unit='$K$', norm=None, xsize=800, lonra=lonra, latra=latra, sub=(4,4,k+4),\
	notext = True, cbar=cbr, margins=(0.005,0.00001,0.005,0.005) )

	# Cal. #
	theta = (90.0-b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)

	hp.projplot(l, b, 'mx', lonlat=True, coord='G', ms=6, mew=2, color='k')

	# for x in pl.frange(l-offset, l+offset, dd):
	# 	for y in pl.frange(b-offset, b+offset, dd):
	# 		cosb = np.cos(b*deg2rad)
	# 		cosy = np.cos(y*deg2rad)
			# if ( (((x-l)**2 + (y-b)**2) <= offset**2) ):
			# if ( (((x-l)**2 + (y-b)**2) <= offset**2) and (((x-l)**2 + (y-b)**2) >= (offset-offset*0.05)**2) ):
				# hp.projtext(x, y, '.', lonlat=True, coord='G')
				# hp.projplot(x, y, 'mx', lonlat=True, coord='G', ms=1, mew=1)

	k = k + 4

mpl.rcParams.update({'font.size':11})
# plt.tight_layout()
# plt.savefig('src_examples.eps', format='eps', dpi=100)
plt.show()