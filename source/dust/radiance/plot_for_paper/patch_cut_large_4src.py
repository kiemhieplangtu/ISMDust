import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class
import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import module            as md

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
def read_nhi_94src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt'):
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
def read_23src_with_OH(fname = '../../oh/result/23src_with_oh.txt'):
	cols = ['src','l','b']
	fmt  = ['s',  'f','f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()

	src  = dat['src']
	xl   = dat['l']
	xb   = dat['b']

	for i in range(len(xl)):
		if (xl[i]>180):
			xl[i] = xl[i]-360.

	return src, xl, xb
##================= MAIN ========================##
pth      = os.getenv("HOME")+'/hdata/dust/'
copth    = os.getenv("HOME")+'/hdata/co/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
sfd_map  = pth + 'lambda_sfd_ebv.fits'
wco_file = copth + 'lambda_wco_dht2001.fits'
co_file  = copth + 'HFI_CompMap_CO-Type2_2048_R1.10.fits'

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
# wco_map   = hp.read_map(wco_file, field = 0)
# wco_map   = hp.read_map(co_file, field = 0)
Xc,Yc,\
wco_map   = md.Wco_map()


dg = 1.5
k  = 0
sc = ['3C18', '3C105', '3C109']
tt = [r'$\tau_{353}$' + '\n'+ '$[10^{-6}]$', 'E(B-V)'+ '\n'+ '$[10^{0}]mag$', 'Radiance'+ '\n'+ '$[10^{-8}]Wm^{-2}sr^{-1}$', 'Temperature'+ '\n'+ '$K$']
tt = [r'$\tau_{353}$', 'E(B-V)', 'Radiance', 'Temperature']

ohsrc, xl, xb = read_23src_with_OH(fname = '../../../oh/result/4src_noCO_noOH_for_plotting.txt')

tmin  = 2.0;  tmax = 8.
tick1 = np.arange(tmin, tmax, 1.)

emin  = 0.03; emax = 0.2
tick2 = np.arange(emin, emax, 0.02)

rmin  = 7.0;  rmax = 42.
tick3 = np.arange(rmin, rmax, 5.)

Tmin  = 17.; Tmax = 22.
tick4 = np.arange(Tmin, Tmax, 0.5)

cbr  = True
xp   = {}
xlra = np.empty([4,2])
xbra = np.empty([4,2])
for i in range(4):
	lonra1 = xl[i]-dg
	lonra2 = xl[i]+dg
	if(lonra1<-180.):
		lonra1 = -180.
		lonra2 = -177.
	if(lonra2>180.):
		lonra2 = 180.
		lonra1 = 177.

	lonra           = [ lonra1, lonra2 ]
	latra           = [ xb[i]-dg, xb[i]+dg ]

	xlra[i][:]      = lonra
	xbra[i][:]      = latra
	xp[ohsrc[i]]    = {}
	xp[ohsrc[i]][0] = hp.cartview(tau_map*1e6, title='tau', coord='G', unit='$(10^{-6})$', norm=None, lonra=lonra, latra=latra, xsize=800,return_projected_map=True )
	xp[ohsrc[i]][1] = hp.cartview(ebv_map*0.86, title='E(B-V)', coord='G', unit='$mag$', norm=None, lonra=lonra, latra=latra, xsize=800, return_projected_map=True )
	xp[ohsrc[i]][2] = hp.cartview(rad_map*1e8,  title='Rad', coord='G', unit='$(10^{-8})Wm^{-2}sr^{-1}$', norm=None, lonra=lonra, latra=latra, xsize=800, return_projected_map=True )
	xp[ohsrc[i]][3] = hp.cartview(tmp_map, title='Temperature', coord='G', unit='$K$', norm=None, lonra=lonra, latra=latra, xsize=800, return_projected_map=True )
	# xp[ohsrc[i]][4] = hp.cartview(wco_map, title='', coord='G', unit='', norm=None, lonra=lonra, latra=latra, xsize=800, return_projected_map=True )


plt.cla()
plt.clf()
plt.close('all')

fig, axarr = plt.subplots(nrows=4, ncols=4,figsize=(15.5,14))
sfmt       = mpl.ticker.ScalarFormatter(useMathText=True) 
# sfmt.set_powerlimits((0, 0))
sfmt.set_scientific(True)
sfmt.set_powerlimits((-1, 1))

cbsize = '1.6%'
fs     = 22
lbsize = 9.5
cbfs   = 18
xysize = 12
cbpad  = 0.02
xmin   = 0.45*3. ## 3-sigma Noise

# for i in range(4):
for i in range(4):
	# xlist = np.linspace(xlra[i][1], xlra[i][0], xp[ohsrc[i]][4].shape[1])    # Create 1-D arrays for x,y dimensions
	# ylist = np.linspace(xbra[i][0], xbra[i][1], xp[ohsrc[i]][4].shape[0])
	# X,Y   = np.meshgrid(xlist, ylist)        # Create 2-D grid xlist,ylist values

	img1  = axarr[0,i].imshow( xp[ohsrc[i]][0], extent=(xlra[i][1],xlra[i][0],xbra[i][0],xbra[i][1]), origin='lower', cmap='jet' )
	axarr[0,i].plot(xl[i], xb[i], '-kx', markersize=7, markeredgewidth=2)

	major_xticks = np.arange(37.5, 40., 1.)
	minor_xticks = np.arange(37.5, 40., 0.5)
	major_yticks = np.arange(59., 62., 1.)
	minor_yticks = np.arange(59., 62., 0.5)
	axarr[i,0].set_xticks(major_xticks)
	axarr[i,0].set_xticks(minor_xticks, minor=True)
	axarr[i,0].set_yticks(major_yticks)
	axarr[i,0].set_yticks(minor_yticks, minor=True)
	axarr[i,0].tick_params(axis='x', labelsize=xysize)
	axarr[i,0].tick_params(axis='y', labelsize=xysize)
	axarr[i,0].tick_params(which='both', width=1.5)
	axarr[i,0].tick_params(which='major', length=6)
	axarr[i,0].tick_params(which='minor', length=3)

	axarr[0,i].set_title(ohsrc[i])
	axarr[0,i].set_xlim(xlra[i][1],xlra[i][0])
	axarr[0,i].set_ylim(xbra[i][0],xbra[i][1])
	axarr[0,0].set_ylabel(r'$\tau_{353}$', fontsize=fs)
	
	divider = make_axes_locatable(axarr[0,i])
	cax1    = divider.append_axes("right", size=cbsize, pad=cbpad)
	cb      =  plt.colorbar(img1, cax=cax1, orientation='vertical')
	if(i==3):
		cb.set_label(label='(10$^{-6}$)',weight='normal',fontsize=cbfs)
		cb.ax.tick_params(labelsize=lbsize)
	else:
		cb.set_label(label='',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)


	

	### Img 2 ###
	img2 = axarr[1,i].imshow( xp[ohsrc[i]][1], extent=(xlra[i][1],xlra[i][0],xbra[i][0],xbra[i][1]), origin='lower' )
	axarr[1,i].plot(xl[i], xb[i], '-kx', markersize=7, markeredgewidth=2)
	
	major_xticks = np.arange(38., 41., 1.)
	minor_xticks = np.arange(38., 41., 0.5)
	major_yticks = np.arange(57., 60., 1.)
	minor_yticks = np.arange(57., 60., 0.5)
	axarr[i,1].set_xticks(major_xticks)
	axarr[i,1].set_xticks(minor_xticks, minor=True)
	axarr[i,1].set_yticks(major_yticks)
	axarr[i,1].set_yticks(minor_yticks, minor=True)
	axarr[i,1].tick_params(axis='x', labelsize=xysize)
	axarr[i,1].tick_params(axis='y', labelsize=xysize)
	axarr[i,1].tick_params(which='both', width=1.5)
	axarr[i,1].tick_params(which='major', length=6)
	axarr[i,1].tick_params(which='minor', length=3)

	axarr[1,i].set_xlim(xlra[i][1],xlra[i][0])
	axarr[1,i].set_ylim(xbra[i][0],xbra[i][1])
	axarr[1,i].tick_params(axis='both', labelsize=xysize)
	axarr[1,0].set_ylabel('E(B-V)', fontsize=fs)

	divider = make_axes_locatable(axarr[1,i])
	cax2    = divider.append_axes("right", size=cbsize, pad=cbpad)
	cb      = plt.colorbar(img2, cax=cax2, orientation='vertical',format=sfmt)
	if(i==3):
		cb.set_label(label='(mag)',weight='normal',fontsize=cbfs)
		cb.ax.tick_params(labelsize=lbsize)
	else:
		cb.set_label(label='',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)
	text = cb.ax.yaxis.get_offset_text()
	text.set_fontsize(lbsize)
	text.set_color("black")
	text.set_position((7.5, 1))

	


	img3 = axarr[2,i].imshow( xp[ohsrc[i]][2], extent=(xlra[i][1],xlra[i][0],xbra[i][0],xbra[i][1]), origin='lower' )
	axarr[2,i].plot(xl[i], xb[i], '-kx', markersize=7, markeredgewidth=2)

	major_xticks = np.arange(28.5, 32., 1.)
	minor_xticks = np.arange(28.5, 32., 0.5)
	major_yticks = np.arange(53.5, 56., 1.)
	minor_yticks = np.arange(53.5, 56., 0.5)
	axarr[i,2].set_xticks(major_xticks)
	axarr[i,2].set_xticks(minor_xticks, minor=True)
	axarr[i,2].set_yticks(major_yticks)
	axarr[i,2].set_yticks(minor_yticks, minor=True)
	axarr[i,2].tick_params(axis='x', labelsize=xysize)
	axarr[i,2].tick_params(axis='y', labelsize=xysize)
	axarr[i,2].tick_params(which='both', width=1.5)
	axarr[i,2].tick_params(which='major', length=6)
	axarr[i,2].tick_params(which='minor', length=3)

	axarr[2,i].set_xlim(xlra[i][1],xlra[i][0])
	axarr[2,i].set_ylim(xbra[i][0],xbra[i][1])
	axarr[2,i].tick_params(axis='both', labelsize=xysize)
	axarr[2,0].set_ylabel('Radiance', fontsize=fs)

	divider = make_axes_locatable(axarr[2,i])
	cax3    = divider.append_axes("right", size=cbsize, pad=cbpad)
	cb      = plt.colorbar(img3, cax=cax3, orientation='vertical')
	if(i==3):
		cb.set_label(label='(10$^{-8}Wm^{-2}sr^{-1}$)',weight='normal',fontsize=cbfs)
		cb.ax.tick_params(labelsize=lbsize)
	else:
		cb.set_label(label='',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)

	### Contour ###
	print xlra[i][0],xlra[i][1],xbra[i][0],xbra[i][1]
	xd,yd,zd = md.patch_Wco_map(Xc,Yc,wco_map, xlra[i][0],xlra[i][1],xbra[i][0],xbra[i][1])
	if(np.all(zd.min() != zd.max())):
		print 'Contour!!!'
		xmax = zd.max()
		xmin = zd.min()
		axarr[2,i].contour(xd, yd, zd, np.arange(xmin, xmax, (xmax-xmin)/4. ), colors='k', linestyles='solid')
		# axarr[1,i].contour(xd, yd, zd, 5, colors='k', linestyles='solid')
	### Contours ###

	


	img4 = axarr[3,i].imshow( xp[ohsrc[i]][3], extent=(xlra[i][1],xlra[i][0],xbra[i][0],xbra[i][1]), origin='lower' )
	axarr[3,i].plot(xl[i], xb[i], '-kx', markersize=7, markeredgewidth=2)
	
	major_xticks = np.arange(36., 38.5, 1.)
	minor_xticks = np.arange(36., 38.5, 0.5)
	major_yticks = np.arange(41.5, 45., 1.)
	minor_yticks = np.arange(41.5, 45., 0.5)
	axarr[i,3].set_xticks(major_xticks)
	axarr[i,3].set_xticks(minor_xticks, minor=True)
	axarr[i,3].set_yticks(major_yticks)
	axarr[i,3].set_yticks(minor_yticks, minor=True)
	axarr[i,3].tick_params(axis='x', labelsize=xysize)
	axarr[i,3].tick_params(axis='y', labelsize=xysize)
	axarr[i,3].tick_params(which='both', width=1.5)
	axarr[i,3].tick_params(which='major', length=6)
	axarr[i,3].tick_params(which='minor', length=3)

	axarr[3,i].set_xlim(xlra[i][1],xlra[i][0])
	axarr[3,i].set_ylim(xbra[i][0],xbra[i][1])
	axarr[3,i].tick_params(axis='both', labelsize=12)
	axarr[3,0].set_ylabel('Temperature', fontsize=fs)

	divider = make_axes_locatable(axarr[3,i])
	cax4    = divider.append_axes("right", size=cbsize, pad=cbpad)
	cb      = plt.colorbar(img4, cax=cax4, orientation='vertical')
	if(i==3):
		cb.set_label(label='(K)',weight='normal',fontsize=cbfs)
		cb.ax.tick_params(labelsize=lbsize)
	else:
		cb.set_label(label='',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)


plt.subplots_adjust(wspace=0.04, hspace=0.01)
fig.tight_layout()
fig.savefig('src_examples.eps', format='eps', bbox_inches='tight', pad_inches=0.05, dpi=100)
fig.show()