import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator
from   restore           import restore

# Find 26 sources with no CO #
#
# params string fname Filename
#
# return void
# 
# Author Van Hiep
##
def read_info_no_co_79(fname = '79src_info.txt'):
	info = read_info()
	#print map(operator.sub, info['b'], info['bc'])

	j=0
	for i in range(0,79):
		if (info['yn'][i] == 0) :
			print('{0}\t{1}\t\t{2}\t\t{3}\t\t{4}\t\t{5}\t\t{6}\t\t{6}'.format(j, info['src'][i], info['l'][i], info['b'][i], info['ra_icrs'][i], info['de_icrs'][i], info['ra_j'][i], info['de_j'][i]))
			j = j + 1


## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict infocd 
 # 
 # version 1/2017
 # Author Van Hiep ##
def read_info_no_co(fname = '../../co12/result/26src_no_co_with_sponge.dat'):
	cols = ['idx','src','l','b','ra_icrs','de_icrs','ra_j','de_j', 'oh', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f','f', 'f',    'f',       's',    's',    'i', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

#================= MAIN ========================#	
deg2rad  = np.pi/180.
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
map_file = pth + 'HFI_CompMap_DustOpacity_2048_R1.10.fits'
# map_file = pth + 'lambda_sfd_ebv.fits'
# map_file = pth + 'lambda_green_dust_map_2d.fits'

info     = read_info_no_co('../../co12/result/26src_no_co_with_sponge.dat')
dbeam    = 3.5/120.0 # Beam = 3.5' -> dbeam = beam/60/2
#dbeam = 0.1

fukui_cf  = 2.10 #2.10e26
planck_cf = 1.18 #1.18e26

## E(B-V) map ##
ci_map = hp.read_map(map_file, verbose=False, field = 2)
nside  = hp.get_nside(ci_map)
res    = hp.nside2resol(nside, arcmin = False)
dd     = res/deg2rad/2.0

#============= For Mask ========= #
offset = 2.0 #degree
ci_map = np.log10(ci_map)
print 'max', ci_map.max()
print 'min', ci_map.min()

#====== For Plotting ======#
title = r'$E(B-V)$'
hp.mollview(ci_map, title='', coord='G', unit=r'$log_{10}(E(B-V))$', norm=None, min=-2.9,max=2.)

for i in range(0,26):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = ci_map[pix]
	print src, l,b,pix, val
	hp.projplot(l, b, 'bo', lonlat=True, coord='G')
	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
	if (b<60):
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, weight='bold')
	else:
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='k', weight='bold')

	if(src == '3C109'):
		hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='b', weight='bold')

mpl.rcParams.update({'font.size':30})
hp.graticule()
plt.grid()
plt.show()

l  = 51.641648
b  = -9.6750019

# Plot cartview a/o mollview #
ll = l
if (l>180):
	ll = ll-360.

offset = 1.
lonr = [51.25, 52.0]
latr = [-10., -9.35]
hp.cartview(ci_map, title='XXX', coord='G', unit='mag', min=0.1,max=0.4,
		norm=None, xsize=800, lonra=lonr, latra=latr, #lonra=[ll-offset,ll+offset], latra=[b-offset,b+offset], min=0, max=0.4,
		return_projected_map=True)

# hp.mollview(tau_map, title=info['src'][i]+'('+str(info['l'][i])+','+str(info['b'][i])+') - '+map_file,
# 	coord='G', unit='', rot=[l,b,0], norm='hist', min=1e-7,max=1e-3, xsize=800)

# Cal. #
hp.projplot(ll, b, 'bo', lonlat=True, coord='G')
hp.projtext(ll, b, ' (' + str(round(ll,2)) + ',' + str(round(b,2)) + ')', lonlat=True, coord='G', fontsize=18, weight='bold')

# theta = (90.0-b)*deg2rad
# phi   = l*deg2rad
# pix   = hp.ang2pix(nside, theta, phi, nest=False)

# for x in pl.frange(l-offset, l+offset, dd):
# 	for y in pl.frange(b-offset, b+offset, dd):
# 		cosb = np.cos(b*deg2rad)
# 		cosy = np.cos(y*deg2rad)
# 		if ( (((x-l)**2 + (y-b)**2) <= offset**2) ):
# 			# hp.projtext(x, y, '.', lonlat=True, coord='G')
# 			hp.projplot(x, y, 'kx', lonlat=True, coord='G')

mpl.rcParams.update({'font.size':30})
plt.show()