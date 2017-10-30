import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import pylab             as pl
import module            as md

from restore             import restore

## Read peaks Info ##
 #
 # params string fname File-name
 # return list guessp Guess-Parameters 
 #
 # version 10/2016 
 # author Nguyen Van Hiep ##
def all_oh_src(fname='oh_yn_src.txt'):
	cols = ['idx','src','l','b','ohyn']
	fmt  = ['i',  's',  'f','f', 'i'  ]
	dat  = restore(fname, 2, cols, fmt)
	info = dat.read(asarray=False)

	return info

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

####============== MAIN ==============####
## low NHI ##
lownhi = md.read_16src_lownhi(fname = '../../oh/result/16src_lowNHI.txt')
lowsrc = lownhi['src']
print ''
print '=> Low NHI src:'
print len(lowsrc)
print lowsrc


## 110 CO Sources  ##
cosrc  = read_co_src(fname='../../co12/data/newly_observed_CO_src.txt')
print ''
print '=> CO src:'
print len(cosrc)
print cosrc

## 23 no CO ##
noCO   = md.read_basic_info_no_co(fname = '../../co12/result/23src_no_co_basic.dat')
noCOsc = noCO['src']
print ''
print '=> NO CO src:'
print len(noCOsc)
print noCOsc

## 16 low-NHI Sources  ##
lownhi = md.read_16src_lownhi(fname = '../result/16src_lowNHI.txt')
lowsrc = lownhi['src']
print ''
print '=> all low NHI src:'
print len(lowsrc)
print lowsrc

## 21SPONGE Claire 30 Sources  ##
sp30sc  = md.read_info_sponge_30src(fname = '../../oh/sponge/sub_data/30src_claire.txt')
spsc    = sp30sc['src']
print ''
print '=> 21-SPONGE src:'
print len(spsc)
print spsc



## OH Sources  ##
ohinfo = all_oh_src('oh_yn_src.txt')
ohsrc  = ohinfo['src']
print ''
print '=> all observed OH sources in datafile:'
print len(ohsrc)
print ohsrc



ohyn   = ohinfo['ohyn']
ohyn   = np.array(ohyn)
fltr   = (ohyn==1)
ohsc   = np.array(ohsrc)
print ''
print '=> all OH-detected sources:'
print sum(1 for x in ohyn if x == 1)
yesOH  = ohsc[ohyn==1]
print yesOH

print ''
print '=> all None-OH-detected sources:'
noOH   = ohsc[ohyn==0]
print sum(1 for x in ohyn if x == 0)
print len(noOH)
print noOH


hi     = md.read_nhi_93src(fname = '../../hi/result/nhi_thin_cnm_wnm_93src_scaled.txt')
hisrc  = hi['src']

for sc in hisrc:
	if(sc in lowsrc):
		print '1'
	else:
		print '-1'

print ''
print '=> all 93 HI sources in (MS+SPONGE):'
print hisrc
print len(hisrc)

msinfo = md.read_ms_79src(fname = '../../hi/data/79_src_nhi_iopscience.txt')
mssrc  = msinfo['src']
print ''
print '=> all 79 HI sources in MS:'
print mssrc
print len(mssrc)


commList1 = md.common_elements(ohsrc,mssrc)
print ''
print 'Common MS and OH sources:'
print len(commList1)


commList2 = md.common_elements(ohsrc,hisrc)
print ''
print 'Common (MS+SP) OH sources:'
print len(commList2)

commList2a = md.common_elements(ohsrc,mssrc)
print ''
print 'MS + OH:'
print len(commList2a)

commList2b = md.common_elements(commList2a,cosrc)
print ''
print 'MS + OH + CO sources:'
print len(commList2b)

commList3 = md.common_elements(commList2,cosrc)
print ''
print 'Common (MS+SP) OH sources and CO sources:'
print len(commList3)


commList4 = md.common_elements(mssrc,cosrc)
print ''
print 'Common MS sources and CO sources:'
print len(commList4)

commList5 = md.common_elements(spsc, ohsrc)
print ''
print 'Common SPONGE sources and OH sources:'
print len(commList5)
print commList5


xl1 = list( set(ohsrc + cosrc) )
print ''
print 'CO + OH sources:'
print len(xl1)
print xl1

xl2 = md.diff_elements(hisrc, xl1)
print ''
print 'Sightlines with HI observation only:'
print len(xl2)
print xl2







print ''
print 'Merged sources OH + HI:'
mergedList = list( set(ohsrc + hisrc) )
print len(mergedList)



cNoOH = md.common_elements(noOH,lowsrc)
print ''
print 'Low NHI and Non-OH-detected sources:'
print len(cNoOH)


print ''
print 'Merged sources OH + CO:'
mergedList = list( set(ohsrc + cosrc) )
print len(mergedList)
for sc in mergedList:
	if(sc in ohsrc):
		if(sc in yesOH):
			print '1'
		else:
			print '0'

	else:
		print '-1'

print ''
print ''
print ''
for sc in mergedList:
	if(sc in lowsrc):
		print '1'
	else:
		print '-1'

























sys.exit()
for sc in mergedList:
	cat = 1
	if( (sc in ohsrc) and (sc in mssrc) ):
		cat = 0

	if( (sc in ohsrc) and (sc not in mssrc) ):
		cat = 2

	print sc, '\t', cat