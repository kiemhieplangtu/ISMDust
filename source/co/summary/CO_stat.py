import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

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
def read_co_src(fname='../../co12/data/newly_observed_CO_src.txt'):
	cols = ['idx', 'src']
	fmt  = ['i', 's']
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read(asarray=False)

	return info['src']

## Read CO src infor ##
 #
 # params string fname File-name
 # return list List of src
 #
 # version 10/2016 
 # author Nguyen Van Hiep ##
def read_info_co_src(fname='../../co12/data/newly_observed_CO_src.txt'):
	cols = ['src', 'coyn']
	fmt  = ['s', 'i']
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read(asarray=False)

	return info['src']

### =============== MAIN ================= ###
observed_src = read_co_src(fname='../data/newly_observed_CO_src.txt')
CO_src       = read_info_co_src(fname='../data/CO_detected_src.dat')
no_CO_src    = read_info_co_src(fname='../data/CO_non_detected_src.dat')

hi           = md.read_nhi_93src(fname = '../../hi/result/nhi_thin_cnm_wnm_93src_scaled.txt')
hisrc        = hi['src']

print ''
print '1) Observed CO sightlines: ', len(observed_src)
print observed_src

print ''
print '2) CO-detected sightlines: ', len(CO_src)
print CO_src

print ''
print '3) CO-non-detected sightlines: ', len(no_CO_src)
print no_CO_src

print ''
print '4) Observed HI sightlines: ', len(hisrc)
hisrc = [x.upper() for x in hisrc]
print hisrc

commList1 = md.common_elements(observed_src,hisrc)
print ''
print 'Common (MS+SP) and CO sources:'
print len(commList1)


diffList1 = md.diff_elements(hisrc, observed_src)
print ''
print 'Diff src from (MS+SP) and CO:'
print len(diffList1)
print diffList1


print ''
for sc in hisrc:
	if(sc in observed_src):
		if(sc in CO_src):
			print '1'
		else:
			print '0'
	else:
		print '-1'


print ''
print ''
print ''
for sc in observed_src:
	if(sc in CO_src):
		print sc, '\t', '1'
	elif(sc in no_CO_src):
		print sc, '\t','0'
	else:
		print sc, '\t','-'

print ''
print ''
print ''
for sc in observed_src:
	if(sc in hisrc):
		print '1'
	else:
		print '0'