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

## Infor for 78 MS sources
src78   = md.read_coord_ms_78src(fname = '../data/78src_radec_lb.txt')
xl      = src78['l']
xb      = src78['b']
src     = src78['src']

msSpec  = md.read_hi_specs(fname = '../data/nhi_opac_specs.txt')

dirs    = os.listdir( '../data/LAB_specs/' )

for i in range(78):
	sc  = src[i]
	zl  = xl[i]
	zb  = xb[i]

	fle = sc+'.txt'
	if(fle not in dirs):
		continue

	vMS    = msSpec[sc]['v']
	texpMS = msSpec[sc]['Texp']

	vLAB,\
	tbLAB  = md.read_LAB_spec(fname = '../data/LAB_specs/'+sc+'.txt')

	# Plot
	fig    = plt.figure(figsize=(12,12))
	ax     = fig.add_subplot(111); #ax.set_rasterized(True)

	plt.plot(vMS, texpMS, 'k-', label='MS')
	plt.plot(vLAB, tbLAB, 'r-', label='LAB')

	plt.title(sc, fontsize = 35)
	plt.ylabel('$T_{b} (K)$', fontsize = 35)
	plt.xlabel('VLSR (km/s)', fontsize = 35)
	plt.tick_params(axis='x', labelsize=20)
	plt.tick_params(axis='y', labelsize=20)

	plt.xlim(-100, 100)

	plt.legend(loc='upper left', fontsize=18)

	plt.savefig("figures/"+sc+'.eps', bbox_inches='tight', pad_inches=0.01, format='eps', dpi=60)
	plt.close(fig)

	# plt.show()