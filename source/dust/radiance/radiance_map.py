import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import module            as md
import operator
from   restore           import restore


#================= MAIN ========================#	
deg2rad  = np.pi/180.
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
# map_file = pth + 'HFI_CompMap_DustOpacity_2048_R1.10.fits'
# map_file = pth + 'lambda_sfd_ebv.fits'
# map_file = pth + 'lambda_green_dust_map_2d.fits'

info     = md.read_19src_noco_nooh(fname = '../../oh/result/19src_noCO_noOH.txt')
lownhi   = md.read_23rc_lownhi(fname = '../../oh/result/16src_lowNHI.txt')
dbeam    = 3.5/120.0 # Beam = 3.5' -> dbeam = beam/60/2
#dbeam = 0.1

# TTYPE1  = 'TAU353  '           / Optical depth at 353GHz                        
# TTYPE2  = 'ERR_TAU '           / Error on optical depth                         
# TTYPE3  = 'EBV     '           / E(B-V) color excess                            
# TTYPE4  = 'RADIANCE'           / Integrated emission                            
# TTYPE5  = 'TEMP    '           / Dust equilibrium temperature                   
# TTYPE6  = 'ERR_TEMP'           / Error on T                                     
# TTYPE7  = 'BETA    '           / Dust emission spectral index                   
# TTYPE8  = 'ERR_BETA'           / error on Beta 

fukui_cf  = 2.10 #2.10e26
planck_cf = 1.18 #1.18e26

## Radiance map ##
r_map  = hp.read_map(map_file, verbose=False, field = 3)
nside  = hp.get_nside(r_map)
res    = hp.nside2resol(nside, arcmin = False)
dd     = res/deg2rad/2.0

#============= For Mask ========= #
offset = 2.0 #degree

#====== For Plotting ======#
hp.cartview(r_map, title='Radiance', coord='G', unit='$Wm^{-2}sr^{-1}$', norm='hist') #, min=0.,max=452)

for i in range(0, len(info['src'])):
	src   = info['src'][i]
	l     = info['l'][i]
	b     = info['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = r_map[pix]
	print src, l,b,pix, val
	hp.projplot(l, b, 'b+', lonlat=True, coord='G')
	# hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
	hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=16, weight='bold')
	# if (b<60):
	# 	hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, weight='bold')
	# else:
	# 	hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='r', weight='bold')

	# if(src == '3C109'):
	# 	hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, color='b', weight='bold')

for i in range(0, len(lownhi['src'])):
	src   = lownhi['src'][i]
	l     = lownhi['l'][i]
	b     = lownhi['b'][i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = r_map[pix]
	print src, l,b,pix, val
	hp.projplot(l, b, 'bo', lonlat=True, coord='G')
	hp.projtext(l, b, src, lonlat=True, coord='G', fontsize=13, weight='bold')

mpl.rcParams.update({'font.size':30})
hp.graticule()
plt.grid()
plt.show()




sys.exit()

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