import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust') # add folder of Class
from common.myImport import *

###================= MAIN ========================###
hiDatPth = const.DATPATH + 'hi/HI4pi/'
map_file = hiDatPth      + 'CAR_-465_465.fits'
hpxfile  = hiDatPth      + 'NHI_HPX.fits'

## Infor for 78 MS sources
src79   = md.read_info_ms_79sc(fname = '../result/nhi_lb_79src_HT03.txt', asarray=True)
xl      = src79['l']
xb      = src79['b']
src     = src79['src']
thin    = src79['thin']
thiner  = src79['thiner']

## N(HI) map ##
HImap   = hp.read_map(hpxfile, field=4)
info    = {'src':src, 'l':xl, 'b':xb}
nhi4pi  = md.get_HI4pi(HImap, info)

# #====== For Plotting ======#
# fig = plt.figure(1, figsize=(32,18))
# hp.mollview(HImap, title='', coord='G', norm='log', sub=(1,1,1), cbar=True, xsize=800, unit='') #min=2.2e-9, max=4.5e-3, format='%0.1e',
# hp.graticule(linestyle=':')

# plt.show()


# sys.exit()

hdulist = fits.open(map_file)

print hdulist.info()
print hdulist[0].header

data    = hdulist[0].data
print('Data shape', data.shape )

hdulist.close()


ddeg    = 0.08333333330000001
xlStart = 180.0 + ddeg
xlEnd   = -xlStart
xbStart = -90.0 - ddeg
xbEnd   = -xbStart


nglBin = 4323
ngbBin = 2163

xglong = np.arange(xlStart, xlEnd, -ddeg)
xglat  = np.arange(xbStart, xbEnd, ddeg)

print xglong.shape
print xglat.shape

# X, Y   = np.meshgrid(xglong,xglat )
# plt.figure(2, figsize=(12,12))
# CS     = plt.contour(X, Y, np.log10(data), 5)
# # CS = plt.contour(np.log10(data), 5)
# plt.clabel(CS, inline=1, fontsize=10)
# plt.title('N*(HI) - HI4pi')

# plt.show()

thin4pi  = np.zeros(79)
for i in range(79):
	sc  = src[i]
	zl  = xl[i]
	zb  = xb[i]

	glId = md.get_index(xglong, zl)
	gbId = md.get_index(xglat, zb)

	# fle = sc+'.txt'
	# if(fle not in dirs):
	# 	continue

	# vMS    = msSpec[sc]['v']
	# texpMS = msSpec[sc]['Texp']

	# line   = open('../data/LAB_specs/'+sc+'.txt', "r").readlines()[2]
	# nhiLAB = line.split(';')[4]
	# nhiLAB = float(nhiLAB)/1e20


	# vLAB,\
	# tbLAB  = md.read_LAB_spec(fname = '../data/LAB_specs/'+sc+'.txt')
	# nhiLab,\
	# er     = cal_optthin_nhi(tbLAB)	

	thin4pi[i] = data[gbId, glId]/1e20

# Plot
fig = plt.figure()
ax  = fig.add_subplot(111); #ax.set_rasterized(True)

plt.plot(thin, nhi4pi, 'ro', label='')
plt.plot(thin, thin4pi, 'k.', label='')
plt.plot([0,50], [0,50], 'k-')

plt.title('HI4pi vs MS', fontsize = 35)
plt.ylabel(r'$N_{HI} - HI4\pi$', fontsize = 35)
plt.xlabel('$N_{HI} - MS$', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)

# plt.xlim(-100, 100)

# plt.legend(loc='upper left', fontsize=18)

# plt.savefig("figures/"+sc+'.eps', bbox_inches='tight', pad_inches=0.01, format='eps', dpi=60)
plt.show()


# Plot for 
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

fig          = plt.figure(figsize=(10,6))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
major_xticks = np.arange(0., 60., 10.)
minor_xticks = np.arange(0., 60., 5.)
major_yticks = np.arange(0., 80., 10.)                                              
minor_yticks = np.arange(0., 80., 5.0)

plt.errorbar(thin, nhi4pi, xerr=thiner*0., yerr=thiner*0., color='r', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, capsize=0, label='From Healpix file')
plt.errorbar(thin, thin4pi, xerr=thiner*0., yerr=thiner*0., color='k', marker='o', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, capsize=0, label='From FITS file')
# plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='None', markerfacecolor='b', markersize=0, label='ODR linear fit')
# plt.fill_between(xfit, mu-sig, mu+sig, color='0.5', alpha=0.5)
plt.plot([0,60], [0,60], 'k--', label='$x = y$')
plt.plot([0,60], [0,60*1.2], 'k--', label='$x = y$')

plt.xlabel(r'$\mathrm{N^*_{HI}\ (HT03)\ [10^{20} cm^{-2}]}$', fontsize=20)
plt.ylabel(r'$\mathrm{N^*_{HI}\ (HI4PI)\ [10^{20} cm^{-2}]}$', fontsize=20)

ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
ax.tick_params(axis='x', labelsize=18, pad=5)
ax.tick_params(axis='y', labelsize=18)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=9)
ax.tick_params(which='minor', length=4)

if(0):
	ax.set_xlim(0.8, 50.)
	ax.set_ylim(0.6, 60.)

	ax.set_xscale('log')
	ax.set_yscale('log')

	plt.savefig('HI4PI_vs_HT03_logscale.eps', pad_inches=0.08, format='eps', dpi=600)

	plt.legend(loc='upper left', fontsize=18)

	plt.tight_layout()
	# plt.savefig('HI4PI_vs_HT03_logscale.eps', bbox_inches='tight', pad_inches=0.08, format='eps', dpi=600)

else:
	ax.set_xlim(-2.0, 50.)
	ax.set_ylim(-2.0, 60.)

	plt.legend(loc='upper left', fontsize=18)

	plt.tight_layout()
	# plt.savefig('HI4PI_vs_HT03.eps', pad_inches=0.08, format='eps', dpi=600)

plt.show()


dirs     = os.listdir( '../data/LAB_specs/' )
thinLAB  = np.zeros(79)
thinLab  = np.zeros(79)
error    = np.zeros(79)
for i in range(79):
	sc  = src[i]
	zl  = xl[i]
	zb  = xb[i]

	line   = open('../data/LAB_specs/'+sc+'.txt', "r").readlines()[2]
	nhiLAB = line.split(';')[4]
	nhiLAB = float(nhiLAB)/1e20


	vLAB,\
	tbLAB  = md.read_LAB_spec(fname = '../data/LAB_specs/'+sc+'.txt')

	thinLAB[i] = nhiLAB


# LAB and HI4PI
fig          = plt.figure(figsize=(10,6))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
major_xticks = np.arange(0., 60., 10.)
minor_xticks = np.arange(0., 60., 5.)
major_yticks = np.arange(0., 80., 10.)                                              
minor_yticks = np.arange(0., 80., 5.0)


plt.errorbar(thinLAB, nhi4pi, xerr=thinLAB*0., yerr=thinLAB*0., color='r', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, capsize=0, label='')
plt.plot([0,60], [0,60], 'k--', label='$x = y$')

plt.xlabel(r'$\mathrm{N^*_{HI}\ (LAB)\ [10^{20} cm^{-2}]}$', fontsize=20)
plt.ylabel(r'$\mathrm{N^*_{HI}\ (HI4PI)\ [10^{20} cm^{-2}]}$', fontsize=20)

ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
ax.tick_params(axis='x', labelsize=18, pad=5)
ax.tick_params(axis='y', labelsize=18)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=9)
ax.tick_params(which='minor', length=4)

plt.savefig('HI4PI_vs_LAB.eps', pad_inches=0.08, format='eps', dpi=600)

plt.show()