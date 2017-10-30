import os, sys
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md

from restore             import restore

ebvInfo = md.read_ebv_av_93src(fname = 'data/ebv_sfd98_sf2011_for_93src.txt', sfd98=False, asarray=True)
ebvSrc  = ebvInfo['src']
ebv     = ebvInfo['ebv']
ebver   = ebvInfo['ebv_er']
ebvsf   = ebvInfo['ebvsf']
ebvsfer = ebvInfo['ebvsf_er']
av      = ebvInfo['av']
zl      = ebvInfo['l']
zb      = ebvInfo['b']

hi      = md.read_nhi_93src(fname = '../../hi/result/nhi_thin_cnm_wnm_93src.txt')
hiSrc   = hi['src']
xl      = hi['l']
xb      = hi['b']

for i in range(len(ebvSrc)):
	# print i, ebvSrc[i], hiSrc[i], 
	# print i, ebvSrc[i], xl[i]==zl[i], xb[i]==zb[i] 
	print i, ebvSrc[i], hiSrc[i], xl[i], xb[i], ebv[i], ebver[i], ebvsf[i], ebvsfer[i], av[i]
