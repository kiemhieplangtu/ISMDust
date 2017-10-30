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
from astropy.io              import fits

###================= MAIN ========================###
deg2rad  = np.pi/180.
pth      = os.getenv("HOME")+'/hdata/hi/HI4pi/'
# map_file = pth + 'CAR_A02.fits'
map_file = pth + 'MOL.fits'
dbeam    = 3.5/120.0 # Beam = 3.5' -> dbeam = beam/60/2


## Radiance map ##
# r_map  = hp.read_map(map_file, verbose=False, field = 0)
# nside  = hp.get_nside(r_map)
# res    = hp.nside2resol(nside, arcmin = False)
# dd     = res/deg2rad/2.0

# #====== For Plotting ======#
# fig = plt.figure(1, figsize=(32,18))
# hp.mollview(r_map, title='', coord='G', norm='log', sub=(1,1,1), cbar=True, xsize=800, min=2.2e-9, max=4.5e-3, format='%0.1e', unit=r'$\tau_{353}$')
# hp.graticule(linestyle=':')

# plt.show()

## Infor for 78 MS sources
src78   = md.read_coord_ms_78src(fname = '../data/78src_radec_lb.txt')
xl      = src78['l']
xb      = src78['b']
src     = src78['src']

msSpec  = md.read_hi_specs(fname = '../data/nhi_opac_specs.txt')



hdulist = fits.open(map_file, mode='denywrite')

print hdulist.info()
header = hdulist[0].header
print header

data     = hdulist[0].data
print('Data shape', data.shape )
print 'Length: '
n = len(data)

hdulist.close()

sys.exit()


print ''
cdelt1 = header['CDELT1']
print cdelt1

glStart = 180.

nglBin = 3891
ngbBin = 1947
vBin   = 945

xglong = np.arange(180., -180., -360./nglBin)
xglat  = np.arange(-90., 90., 180./ngbBin)
vlsr   = np.arange(-470., 470., 940./vBin)

print nglBin

for i in range(78):
	sc  = src[i]
	zl  = xl[i]
	zb  = xb[i]

	xlId = md.get_index(xglong, zl)
	xbId = md.get_index(xglat, zb)

	print '-----------'
	print 'Src: ', sc
	print 'L,B & indices'
	print zl, xlId
	print zb, xbId
	print ''

	T      = data[:, xbId, xlId]

	print T

	vMS    = msSpec[sc]['v']
	texpMS = msSpec[sc]['Texp']

	# Plot
	fig    = plt.figure(figsize=(12,12))
	ax     = fig.add_subplot(111); #ax.set_rasterized(True)

	plt.plot(vMS, texpMS, 'k-', label='MS')
	plt.plot(vlsr, T, 'r-', label='HI4pi')

	plt.title(sc, fontsize = 35)
	plt.ylabel('$T_{b} (K)$', fontsize = 35)
	plt.xlabel('VLSR (km/s)', fontsize = 35)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)

	plt.legend(loc='upper left', fontsize=18)

	plt.show()

