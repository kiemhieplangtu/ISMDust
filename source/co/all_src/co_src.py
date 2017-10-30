import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import module            as md

from restore             import restore

## Read CO src in file1 ##
 #
 # params string fname File-name
 # return list List of src
 #
 # version 10/2016 
 # author Nguyen Van Hiep ##
def read_co_src_in_datafile(fname='../data/psm2016_millennium_11src.txt'):
	cols = ['src']
	fmt  = ['s']
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read(asarray=False)

	return info['src']

####============== MAIN ==============####

src1 = read_co_src_in_datafile(fname='../data/psm2016_millennium_11src.txt')
src2 = read_co_src_in_datafile(fname='../data/pmodlh14_161src.txt')

src  = read_co_src_in_datafile(fname='../data/old_data_sources.txt')
srcN = read_co_src_in_datafile(fname='../data/src_from_ningyuEmail.txt')

print ''
print '=> all observed CO sources in datafile2:'
src = list(set(src))
print src
print len(src)

commonList = md.common_elements(src,srcN)
print ''
print 'Common CO sources:'
print len(commonList)

diffsrc    = md.diff_elements(src,srcN)
print ''
print 'Diff CO sources:'
print len(diffsrc)
print diffsrc


print ''
print '=> all observed CO sources in datafile1:'
print src1
print len(src1)

print ''
print '=> all observed CO sources in datafile2:'
src2 = list(set(src2))
print src2
print len(src2)


commonList = md.common_elements(src1,src2)
print ''
print 'Common CO sources:'
print len(commonList)

print ''
print 'Merged sources:'
mergedList = list( set(src1 + src2) )
print len(mergedList)

# print ''
# print 'Merged sources:'
# for i in range(len(mergedList)):
# 	print i, mergedList[i]

print ''
print 'Sources:'
for i in range(len(srcN)):
	if(srcN[i] in mergedList):
		print i+1, srcN[i], '  - yes'
	else:
		print i+1, srcN[i], '  - no'