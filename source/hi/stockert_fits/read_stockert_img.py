import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class
import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp

import operator
from restore             import restore
from astropy.io import fits

## Read info of 26 no-CO sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 11/2016
 # Author Van Hiep ##
def read_23oh_src(fname = '../oh/result/23src_with_oh.txt'):
	cols = ['src','l','b']
	fmt  = ['s',  'f','f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

#================= MAIN ========================#	
## Define constants ##
deg2rad = np.pi/180.
beam    = 3.5
dbeam   = beam/120.0 # Beam = 3.5' -> dbeam = beam/60/2
dbeam   = 0.5
edge    = dbeam/40. # To plot a rectangle about the source
# info    = read_23oh_src('../oh/result/23src_with_oh.txt')

print 'mmm: ', 2.8+6.32*(1420./1665.402)**2.8
print 'mmm: ', 2.8+22.6*(408./1665.402)**2.8

# hdulist = fits.open('data/old/fitsp042820_17min.bin')
# hdu = hdulist[0]

# print hdu.data.shape
# print hdu.data
# print np.mean(hdu.data)
# # print hdu.header
# plt.imshow(hdu.data[0,:,:], origin='lower')
# plt.show()

hdulist = fits.open('data/fitsp042820.bin')
hdu = hdulist[0]

print hdu.data.shape
print hdu.data
print np.mean(hdu.data)
# print hdu.header
plt.imshow(hdu.data[0,:,:], origin='lower')
plt.show()

# Product Name
# CHIPASS 1.4 GHz Continuum Map
# Coord. System
# Galactic
# Projection Type
# HEALPix, nested, res 10 (Nside=1024)
# Resolution
# 14.4 arcmin
# Original Data Source ATNF

# himap = hp.read_map('fits30587.bin', field = 0, h=False)
# nside = hp.get_nside(himap)

## Color map ##
# cmap  = plt.cm.get_cmap('cool')
# hp.mollview(himap, title='21cm radio Continuum', coord='G', unit='', rot=[0,0,0], norm='log', xsize=800, min=1000.0, cmap=cmap)

# for i in range(0,len(info['src'])):
# 	src   = info['src'][i]
# 	l     = info['l'][i]
# 	b     = info['b'][i]

# 	theta = (90.0 - b)*deg2rad
# 	phi   = l*deg2rad
# 	pix   = hp.ang2pix(nside, theta, phi, nest=False)
# 	val   = himap[pix]
# 	print src, l,b,pix, val
# 	hp.projplot(l, b, 'bx', lonlat=True, coord='G')
# 	hp.projtext(l, b, src+','+str(l)+','+str(b), lonlat=True, coord='G')
# 	# hp.projtext(l, b, src+','+str(val), lonlat=True, coord='G')
# 	# hp.projtext(l, b, src, lonlat=True, coord='G')

# hp.graticule()
# plt.show()