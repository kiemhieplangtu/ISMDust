import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import module            as md

from restore             import restore

## Read 110 CO src ##
 #
 # params string fname File-name
 # return list List of src
 #
 # version 10/2016 
 # author Nguyen Van Hiep ##
def read_src(fname='OH_observed_src.dat'):
	cols = ['src']
	fmt  = ['s']
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read(asarray=False)

	return info['src']

### =============== MAIN ================= ###
observed_src = read_src(fname='OH_observed_src.dat')
OH_src       = read_src(fname='OH_detected_src.dat')
no_OH_src    = read_src(fname='OH_non_detected_src.dat')
hi           = md.read_nhi_93src(fname = '../../hi/result/nhi_thin_cnm_wnm_93src_scaled.txt')
hisrc        = hi['src']

print ''
print '1) Observed sightlines: ', len(observed_src)
print observed_src

print ''
print '2) OH detected sightlines: ', len(OH_src)
print OH_src

print ''
print '3) OH non-detected sightlines: ', len(no_OH_src)
print no_OH_src

print ''
for sc in hisrc:
	if(sc in observed_src):
		if(sc in OH_src):
			print '1'
		else:
			print '0'
	else:
		print '-1'


print ''
print ''
print ''
for sc in observed_src:
	if(sc in OH_src):
		print sc, '\t', 1
	else:
		print sc, '\t', 0

print ''
print ''
print ''
for sc in observed_src:
	if(sc in hisrc):
		print 1
	else:
		print 0