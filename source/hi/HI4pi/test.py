import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust') # add folder of Class

from common.myImport import *

###================= MAIN ========================###
deg2rad  = np.pi/180.
pth      = os.getenv("HOME")+'/hdata/hi/HI4pi/'
# map_file = pth + 'CAR_A02.fits'
map_file = pth + 'HPX_001.fits'
dbeam    = 3.5/120.0 # Beam = 3.5' -> dbeam = beam/60/2


## Radiance map ##
r_map  = hp.read_map(map_file, verbose=False, field = 0)
# nside  = hp.get_nside(r_map)
# res    = hp.nside2resol(nside, arcmin = False)
# dd     = res/deg2rad/2.0

# #====== For Plotting ======#
fig = plt.figure(1, figsize=(32,18))
hp.mollview(r_map, title='', coord='G', norm='log', sub=(1,1,1), cbar=True, xsize=800, min=2.2e-9, max=4.5e-3, format='%0.1e', unit=r'$\tau_{353}$')
hp.graticule(linestyle=':')

plt.show()


hdulist = fits.open(map_file)

print hdulist.info()
print hdulist[0].header

# data     = hdulist[0].data
# print('WCO data shape', data.shape )
# print data

# T = data[:, 50, 50]

# print T

# plt.plot(T)
# plt.show()

hdulist.close()