import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class
import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import module            as md

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
def read_23src_with_OH(fname = '../../oh/result/23src_with_oh_0.txt'):
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
# info  = read_info_no_co('26src_no_co_info.dat')
info = read_nhi_94src(fname = '../../../../hi/result/nhi_thin_cnm_wnm_94src_sponge_prior.txt')
src  = info['src']

## 4 first sources ##
ohsrc, xl, xb = read_23src_with_OH(fname = '../../../../oh/result/23src_with_oh_0.dat')

## Maps ##
pth      = os.getenv("HOME")+'/hdata/dust/'
copth    = os.getenv("HOME")+'/hdata/co/'
map_file = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
sfd_map  = pth + 'lambda_sfd_ebv.fits'
wco_file = copth + 'lambda_wco_dht2001.fits'
co_file  = copth + 'HFI_CompMap_CO-Type2_2048_R1.10.fits'


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

## Width of patche ##
dg    = 1.5

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
	lonra1     = xl[i]-dg
	lonra2     = xl[i]+dg

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
	xp[ohsrc[i]][1] = hp.cartview(ebv_map, title='E(B-V)', coord='G', unit='$mag$', norm=None, lonra=lonra, latra=latra, xsize=800, return_projected_map=True )
	xp[ohsrc[i]][2] = hp.cartview(rad_map*1e8,  title='Rad', coord='G', unit='$(10^{-8})Wm^{-2}sr^{-1}$', norm=None, lonra=lonra, latra=latra, xsize=800, return_projected_map=True )
	xp[ohsrc[i]][3] = hp.cartview(tmp_map, title='Temperature', coord='G', unit='$K$', norm=None, lonra=lonra, latra=latra, xsize=800, return_projected_map=True )
	# xp[ohsrc[i]][4] = hp.cartview(wco_map, title='', coord='G', unit='', norm=None, lonra=lonra, latra=latra, xsize=800, return_projected_map=True )


plt.cla()
plt.clf()
plt.close('all')

sfmt   = mpl.ticker.ScalarFormatter(useMathText=True) 
# sfmt.set_powerlimits((0, 0))
sfmt.set_scientific(True)
sfmt.set_powerlimits((-2, 2))

figsz  = (7,7)
fts    = 32
cbsize = '1.6%'
lbsize = 14
cbpad  = 0.02
xmin   = 0.45*3. ## 3-sigma Noise
dpi    = 50
dp     = 0.02
twdth  = 1.5
mjlen  = 6
mnlen  = 3

for i in range(4):
	## Fig c1 ##
	fig  = plt.figure(1, figsize=figsz)
	ax   = fig.add_subplot(111)
	img1 = plt.imshow( xp[ohsrc[i]][0], extent=(xlra[i][1],xlra[i][0],xbra[i][0],xbra[i][1]), origin='lower', cmap='jet' )
	plt.plot(xl[i], xb[i], '-kx', markersize=7, markeredgewidth=2)

	plt.title(ohsrc[i])
	plt.xlim(xlra[i][1],xlra[i][0])
	plt.ylim(xbra[i][0],xbra[i][1])
	if(i==0):
		plt.ylabel(r'$\tau_{353}$', fontsize=fts)
	plt.tick_params(axis='both', labelsize=lbsize)
	ax.tick_params(which='both', width=twdth)
	ax.tick_params(which='major', length=mjlen)
	ax.tick_params(which='minor', length=mnlen)

	
	divider = make_axes_locatable(ax)
	cax1    = divider.append_axes("right", size=cbsize, pad=cbpad)
	# cb      =  plt.colorbar(img1, cax=cax1, orientation='vertical', ticks=np.arange(np.min(xp[ohsrc[i]][0]), np.max(xp[ohsrc[i]][0]), (np.max(xp[ohsrc[i]][0])-np.min(xp[ohsrc[i]][0]))/5.),format=sfmt)
	cb      =  plt.colorbar(img1, cax=cax1, orientation='vertical')
	if(i==3):
		cb.set_label(label='(10$^{-6}$)',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)
	else:
		cb.set_label(label='',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)
	plt.tight_layout()
	plt.savefig('src_eg_apd0_r0c'+str(i)+'.eps', format='eps', bbox_inches='tight', pad_inches=dp, dpi=dpi)
	plt.cla()
	plt.clf()


	
	## Fig c2 ##
	fig  = plt.figure(2, figsize=figsz)
	ax   = fig.add_subplot(111)
	img2 = ax.imshow( xp[ohsrc[i]][1], extent=(xlra[i][1],xlra[i][0],xbra[i][0],xbra[i][1]), origin='lower' )
	plt.plot(xl[i], xb[i], '-kx', markersize=7, markeredgewidth=2)
	plt.xlim(xlra[i][1],xlra[i][0])
	plt.ylim(xbra[i][0],xbra[i][1])
	plt.tick_params(axis='both', labelsize=lbsize)
	ax.tick_params(which='both', width=twdth)
	ax.tick_params(which='major', length=mjlen)
	ax.tick_params(which='minor', length=mnlen)

	if(i==0):
		plt.ylabel('E(B-V)', fontsize=fts)

	divider = make_axes_locatable(ax)
	cax2    = divider.append_axes("right", size=cbsize, pad=cbpad)
	cb      = plt.colorbar(img2, cax=cax2, orientation='vertical',format=sfmt)
	if(i==3):
		cb.set_label(label='(mag)',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)
	else:
		cb.set_label(label='',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)
	text = cb.ax.yaxis.get_offset_text()
	text.set_fontsize(8)
	text.set_color("black")
	text.set_position((7.5, 1))
	plt.tight_layout()
	plt.savefig('src_eg_apd0_r1c'+str(i)+'.eps', format='eps', bbox_inches='tight', pad_inches=dp, dpi=dpi)
	plt.cla()
	plt.clf()

	

	## Fig c3 ##
	fig  = plt.figure(3, figsize=figsz)
	ax   = fig.add_subplot(111)
	img3 = ax.imshow( xp[ohsrc[i]][2], extent=(xlra[i][1],xlra[i][0],xbra[i][0],xbra[i][1]), origin='lower' )
	plt.plot(xl[i], xb[i], '-kx', markersize=7, markeredgewidth=2)
	plt.xlim(xlra[i][1],xlra[i][0])
	plt.ylim(xbra[i][0],xbra[i][1])
	plt.tick_params(axis='both', labelsize=lbsize)
	ax.tick_params(which='both', width=twdth)
	ax.tick_params(which='major', length=mjlen)
	ax.tick_params(which='minor', length=mnlen)

	if(i==0):
		plt.ylabel('Radiance', fontsize=fts)

	divider = make_axes_locatable(ax)
	cax3    = divider.append_axes("right", size=cbsize, pad=cbpad)
	cb      = plt.colorbar(img3, cax=cax3, orientation='vertical')
	if(i==3):
		cb.set_label(label='(10$^{-8}Wm^{-2}sr^{-1}$)',weight='normal')
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
		ax.contour(xd, yd, zd, np.arange(xmin, xmax, (xmax-xmin)/6. ), colors='k', linestyles='solid')
		# axarr[1,i].contour(xd, yd, zd, 5, colors='k', linestyles='solid')
	### Contours ###
	plt.tight_layout()
	plt.savefig('src_eg_apd0_r2c'+str(i)+'.eps', format='eps', bbox_inches='tight', pad_inches=dp, dpi=dpi)
	plt.cla()
	plt.clf()

	

	## Fig c4 ##
	fig  = plt.figure(4, figsize=figsz)
	ax   = fig.add_subplot(111)
	img4 = ax.imshow( xp[ohsrc[i]][3], extent=(xlra[i][1],xlra[i][0],xbra[i][0],xbra[i][1]), origin='lower' )
	plt.plot(xl[i], xb[i], '-kx', markersize=7, markeredgewidth=2)
	plt.xlim(xlra[i][1],xlra[i][0])
	plt.ylim(xbra[i][0],xbra[i][1])
	plt.tick_params(axis='both', labelsize=lbsize)
	ax.tick_params(which='both', width=twdth)
	ax.tick_params(which='major', length=mjlen)
	ax.tick_params(which='minor', length=mnlen)

	if(i==0):
		plt.ylabel('Temperature', fontsize=fts)

	divider = make_axes_locatable(ax)
	cax4    = divider.append_axes("right", size=cbsize, pad=cbpad)
	cb      = plt.colorbar(img4, cax=cax4, orientation='vertical')
	if(i==3):
		cb.set_label(label='(K)',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)
	else:
		cb.set_label(label='',weight='normal')
		cb.ax.tick_params(labelsize=lbsize)
	plt.tight_layout()
	plt.savefig('src_eg_apd0_r3c'+str(i)+'.eps', format='eps', bbox_inches='tight', pad_inches=dp, dpi=dpi)
	plt.cla()
	plt.clf()


# plt.subplots_adjust(wspace=0.04, hspace=0.01)
# fig.tight_layout()
# fig.savefig('src_examples_apd0.eps', format='eps', bbox_inches='tight', pad_inches=0.05, dpi=100)
# fig.show()