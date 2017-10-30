import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md

from restore             import restore

###================= MAIN ========================####
# Cal tau353 from map or not
readmap = False

## Filename of the maps
pth     = os.getenv("HOME")+'/hdata/dust/'
mapFile = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

tauMap   = hp.read_map(mapFile, field = 0)
tauerMap = hp.read_map(mapFile, field = 1)
rMap     = hp.read_map(mapFile, field = 3)

## Read Infor for 93 src
dat   = md.read_info_93src_csv(fname = '../../../doc/observed_src.csv',asarray=True)
xsc   = dat['src']
xl    = dat['l']
xb    = dat['b']
nhi   = dat['nhi']   ## Already in ascending order
nhier = dat['nhier']
xOK   = dat['ok']
tsg67 = dat['tsig67']
ebv   = dat['ebv']
ebver = dat['ebver']
n     = len(xsc)

# fltr  = np.extract([xOK == 1], xOK)
# xsc   = np.extract([xOK == 1], xsc)
# xl    = np.extract([xOK == 1], xl)
# xb    = np.extract([xOK == 1], xb)
# nhi   = np.extract([xOK == 1], nhi)
# nhier = np.extract([xOK == 1], nhier)
# tsg67 = np.extract([xOK == 1], tsg67)

# to dict. of Infor
infor        = {}
infor['src'] = xsc
infor['l']   = xl
infor['b']   = xb

nh1, nh1er, t353, t353er = md.get_nh_from_tau353(tauMap, tauerMap, infor)
nh2, nh2er, ebv, ebver   = md.get_nh_from_ebv(ebv, ebver)
nh3, nh3er, R, Rer       = md.get_nh_from_radiance(tauMap, tauerMap, rMap, infor)

for i in range(n):
	print i, xsc[i], xl[i], xb[i], nhi[i], nhier[i], t353[i]*1e6, t353er[i]*1e6, nh1[i], nh1er[i], ebv[i], ebver[i], nh2[i], nh2er[i], R[i]*1e11, Rer[i]*1e11, nh3[i], nh3er[i]