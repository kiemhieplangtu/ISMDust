import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class
import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl

import operator
from restore             import restore

## Read info of E(B-V) for OH sources only #
 # src | l | b | E(B-V)SFD_online | err | E(B-V)S&F2011 | err
 #
 # params string fname Filename
 # return dict info
 # 
 # version 04/2017
 # Author Van Hiep ##
def read_ebv_for_oh_src(fname = 'data/ebv_sfd98_sf2011_for_oh_src.txt', sfd98=False):
	cols = ['idx','src','l', 'b', 'ebv','ebv_er', 'ebvsf','ebvsf_er']
	fmt  = ['i',  's',  'f', 'f', 'f',  'f'     , 'f',     'f'      ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()

	if(sfd98):
		return dat['ebv'], dat['ebv_er'], dat['src'], dat['l'], dat['b']
	else:
		return dat['ebvsf'], dat['ebvsf_er'], dat['src'], dat['l'], dat['b']

##================= MAIN ========================##
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
sfd_map  = pth + 'lambda_sfd_ebv.fits'
quantity = 'tau353'

# info  = read_info_no_co('26src_no_co_info.dat')
ebv, ebver, src, xl, xb = read_ebv_for_oh_src(fname = '../dust/ebv2nh/data/ebv_sfd98_sf2011_for_oh_src.txt', sfd98=True)

# Define the width of area #
beam   = 5.0            # Beam = 5' = map_resolution
dbeam  = beam/120.0     # Beam = 3.5' -> dbeam = beam/60/2 in degree
offset = dbeam          # degree

# Define constants #
nside     = 2048
deg2rad   = np.pi/180.

tau_map   = hp.read_map(map_file, field = 0)
ebv_map   = hp.read_map(map_file, field = 2)
rad_map   = hp.read_map(map_file, field = 3)
tmp_map   = hp.read_map(map_file, field = 4)
nside     = hp.get_nside(tau_map)
res       = hp.nside2resol(nside, arcmin=False)
dd        = res/deg2rad/15.0

dg = 1.0
k  = 0
sc = ['3C207', '3C454.3', '3C18']
sc = ['3C410', 'T0526+24', 'T0629+10']
tt = [r'$\tau_{353}$' + '\n'+ '$[10^{-6}]$', 'E(B-V)'+ '\n'+ '$[10^{0}]mag$', 'Radiance'+ '\n'+ '$[10^{-8}]Wm^{-2}sr^{-1}$', 'Temperature'+ '\n'+ '$K$']
tt = [r'$\tau_{353}$', 'E(B-V)', 'Radiance', 'Temperature']
## Plot
fig = plt.figure(figsize=(10,9.))
fig.set_rasterized(True)
for i in range(len(src)):
	if k > 0:
		tt =['','','','']

	# cbr = True
	# if(k<8):
	cbr = True

	if( src[i] not in sc ) :
		continue

	l  = xl[i]
	b  = xb[i]

	print src[i], l, b
	# Plot cartview a/o mollview #
	ll = l
	if (l>180):
		ll = ll-360.

	hp.cartview(tau_map*1e6, title=tt[0], coord='G', unit=src[i], norm=None, xsize=800, lonra=[ll-dg,ll+dg], latra=[b-dg,b+dg], sub=(3,4,k+1),\
	notext = True, cbar=cbr, margins=(0.005,0.00001,0.005,0.005) ) ##  min=2.5, max=65.0, 
	# hp.projtext(ll, b, '  (' + str(round(ll,2)) + ', ' + str(round(b,2)) + ')', lonlat=True, coord='G', fontsize=15, weight='bold', color='y')
	hp.cartview(ebv_map, title=tt[1], coord='G', unit=src[i], norm=None, xsize=800, lonra=[ll-dg,ll+dg], latra=[b-dg,b+dg], sub=(3,4,k+2),\
	notext = True, cbar=cbr, margins=(0.005,0.00001,0.005,0.005))
	hp.cartview(rad_map*1e8, title=tt[2], coord='G', unit=src[i], norm=None, xsize=800, lonra=[ll-dg,ll+dg], latra=[b-dg,b+dg], sub=(3,4,k+3),\
	notext = True, cbar=cbr, margins=(0.005,0.00001,0.005,0.005) )  ## min=8.5, max=160.0, 
	hp.cartview(tmp_map, title=tt[3], coord='G', unit=src[i], norm=None, xsize=800, lonra=[ll-dg,ll+dg], latra=[b-dg,b+dg], sub=(3,4,k+4),\
	notext = True, cbar=cbr, margins=(0.005,0.00001,0.005,0.005) )  ## min=16.5, max=22.5,

	# min=7.03e-6, max=1.24e-5 )
	# hp.projplot(ll, b, 'bo', lonlat=True, coord='G')
	# hp.projtext(ll, b, '  (' + str(round(ll,2)) + ', ' + str(round(b,2)) + ')', lonlat=True, coord='G', fontsize=18, weight='bold', color='y')

	# hp.mollview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file,
	# 	coord='G', unit='', rot=[l,b,0], norm='hist', min=1e-7,max=1e-3, xsize=800)

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