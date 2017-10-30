import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import module            as md
import operator

from   restore               import restore
from mpl_toolkits.axes_grid1 import make_axes_locatable

###================= MAIN ========================###
info = md.read_info_93src_csv(fname = '../../../doc/observed_src.csv',asarray=False)
xsc  = info['src']
xl   = info['l']
xb   = info['b']
ohyn = info['oh']
coyn = info['co']
atom = info['ok']

deg2rad  = np.pi/180.
pth      = os.getenv("HOME")+'/hdata/dust/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
dbeam    = 3.5/120.0 # Beam = 3.5' -> dbeam = beam/60/2

cbsize = '0.005%'
lbsize = 7
cbpad  = 2

## Radiance map ##
r_map  = hp.read_map(map_file, verbose=False, field = 0)
nside  = hp.get_nside(r_map)
res    = hp.nside2resol(nside, arcmin = False)
dd     = res/deg2rad/2.0

#====== For Plotting ======#
fig = plt.figure(1, figsize=(32,18))
hp.mollview(r_map, title='', coord='G', norm='log', sub=(1,1,1), cbar=True, xsize=800, min=2.2e-9, max=4.5e-3, format='%0.1e', unit=r'$\tau_{353}$')
hp.graticule(linestyle=':')

hp.projplot(0., 0., color='k', marker='x', markersize=8, lonlat=True, coord='G', markerfacecolor="None", markeredgecolor='k', markeredgewidth=2)

for i in range(len(xsc)):
	src   = xsc[i]
	l     = xl[i]
	b     = xb[i]

	theta = (90.0 - b)*deg2rad
	phi   = l*deg2rad
	pix   = hp.ang2pix(nside, theta, phi, nest=False)
	val   = r_map[pix]
	print src, l,b,pix, val



	if(ohyn[i] == 1):
		hp.projplot(l, b, color='r', marker='o', markersize=16, markeredgecolor='r', mfc="None", markeredgewidth=2, lonlat=True, coord='G')
	if(ohyn[i] == 0 ):
		hp.projplot(l, b, color='k', marker='o', markersize=16, markeredgecolor='k', mfc="None", markeredgewidth=2, lonlat=True, coord='G')

	if(coyn[i] == 1):
		hp.projplot(l, b, color='r', marker='^', markersize=7, markeredgecolor='r', mfc="r", markeredgewidth=1, lonlat=True, coord='G')
	if(coyn[i] == 0):
		hp.projplot(l, b, color='k', marker='^', markersize=7, markeredgecolor='k', mfc="k", markeredgewidth=1, lonlat=True, coord='G')

	if(atom[i]==1):
		hp.projplot(l, b, color='r', marker='s', markersize=25, markeredgecolor='r', mfc="None", markeredgewidth=2.5, lonlat=True, coord='G')
		# hp.projtext(l, b, ' '+src, lonlat=True, coord='G', fontsize=12, weight='bold', color='k')
	else:
		hp.projplot(l, b, color='k', marker='s', markersize=23, markeredgecolor='k', mfc="None", markeredgewidth=2, lonlat=True, coord='G')

	if(src == "3C132"):
		hp.projtext(l+1.1, b-2.5, ' '+src, lonlat=True, coord='G', fontsize=12, weight='bold', color='k')

# 	theta = (90.0 - b)*deg2rad
# 	phi   = l*deg2rad
# 	pix   = hp.ang2pix(nside, theta, phi, nest=False)
# 	val   = r_map[pix]
# 	print src, l,b,pix, val
# 	hp.projplot(l, b, color='r', marker='^', markersize=10, lonlat=True, coord='G')

mpl.rcParams.update({'font.size':16})
fig.axes[1].texts[0].set_fontsize(32)

plt.savefig('src_locations.eps', bbox_inches='tight', pad_inches=0.01, format='eps', dpi=60)
# plt.show()




sys.exit()