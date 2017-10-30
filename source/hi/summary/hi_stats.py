import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add fnewer of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import module            as md

from restore             import restore

hi93sc = md.read_nhi_93src(fname = '../result/nhi_thin_cnm_wnm_93src.txt')
sc93   = hi93sc['src']
id93   = hi93sc['indx']
xl93   = hi93sc['l']
xb93   = hi93sc['b']
nhi    = hi93sc['nhi']
nhier  = hi93sc['nhi_er']
thin   = hi93sc['thin']
thiner = hi93sc['thin_er']
cnm    = hi93sc['cnm']
cnmer  = hi93sc['cnm_er']
wnm    = hi93sc['wnm']
wnmer  = hi93sc['wnm_er']


## print out NHI infor ###
for i in range(len(sc93)):
	print i+1, sc93[i], xl93[i], xb93[i], nhi[i], nhier[i], thin[i], thiner[i], cnm[i], cnmer[i], wnm[i], wnmer[i]


### Low N(HI), N(HI)<3e20
for i in range(len(sc93)):
	if(nhi[i] < 1.0):
		print 1
	else:
		print 0
