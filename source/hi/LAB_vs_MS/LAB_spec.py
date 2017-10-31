import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

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
pth      = os.getenv("HOME")+'/hdata/hi/LAB/'
map_file = pth + 'lab250.fit'
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

hdulist = fits.open(map_file)

# print hdulist.info()
print hdulist[0].header
data    = hdulist[0].data

glStart = hdulist[0].header['CRVAL1']
gdel    = hdulist[0].header['CDELT1']
nglBin  = hdulist[0].header['NAXIS1']
glEnd   = glStart + nglBin*gdel

print 'Starting G.Longitude:'
print glStart
print 'Delta GLong: '
print gdel
print 'N bins for GLong:'
print nglBin

bStart  = hdulist[0].header['CRVAL2']
bdel    = hdulist[0].header['CDELT2']
nbBin   = hdulist[0].header['NAXIS2']
bEnd    = bStart + nbBin*bdel

print ''
print 'Starting G.Latitude:'
print bStart
print 'Delta GLat: '
print bdel
print 'N bins for GLat:'
print nbBin

nvBin  = hdulist[0].header['NAXIS3']
vStart = hdulist[0].header['CRVAL3']  # m/s
vdel   = hdulist[0].header['CDELT3']  # m/s
vEnd   = vStart + nvBin*vdel          # m/s

hdulist.close()

print ''
print 'Starting Vel:'
print vStart
print 'Delta Vel: '
print vdel
print 'N bins for Vel:'
print nvBin

xglong = np.arange(glStart, glEnd, gdel)
xglat  = np.arange(bStart, bEnd, bdel)
vlsr   = np.arange(vStart, vEnd, vdel)
vlsr   = vlsr/1e3

for i in range(78):
	sc  = src[i]
	zl  = xl[i]
	zb  = xb[i]

	xlId = md.get_index(xglong, zl)
	xbId = md.get_index(xglat, zb)

	print 'L,B indices'
	print xlId
	print xbId

	T      = data[:, xbId, xlId]

	vMS    = msSpec[sc]['v']
	texpMS = msSpec[sc]['Texp']

	# Plot
	fig    = plt.figure(figsize=(12,12))
	ax     = fig.add_subplot(111); #ax.set_rasterized(True)

	plt.plot(vMS, texpMS, 'k-', label='MS')
	plt.plot(vlsr, T, 'r-', label='LAB')

	plt.title(sc, fontsize = 35)
	plt.ylabel('$T_{b} (K)$', fontsize = 35)
	plt.xlabel('VLSR (km/s)', fontsize = 35)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)

	plt.legend(loc='upper left', fontsize=18)

	plt.show()

read_LAB_spec(fname = '../source/hi/data//LAB_specs/abs.txt')