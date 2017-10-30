import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import module            as md
import operator
from   restore           import restore

## Read info of each source to dictionary #
 #
 # params dict dat Data
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read2dict(dat):
	sc     = dat['src']
	xl     = dat['l']
	xb     = dat['b']
	hi     = dat['nhi']
	hier   = dat['nhi_er']
	thin   = dat['thin']
	thiner = dat['thin_er']
	cnm    = dat['cnm']
	cnmer  = dat['cnm_er']
	wnm    = dat['wnm']
	wnmer  = dat['wnm_er']

	ret = {}
	for i in range(len(sc)):
		ret[sc[i]]            = {}
		ret[sc[i]]['l']       = xl[i]
		ret[sc[i]]['b']       = xb[i]
		ret[sc[i]]['nhi']     = hi[i]
		ret[sc[i]]['nhi_er']  = hier[i]
		ret[sc[i]]['thin']    = thin[i]
		ret[sc[i]]['thin_er'] = thiner[i]
		ret[sc[i]]['cnm']     = cnm[i]
		ret[sc[i]]['cnm_er']  = cnmer[i]
		ret[sc[i]]['wnm']     = wnm[i]
		ret[sc[i]]['wnm_er']  = wnmer[i]

	return ret

## ======== MAIN ============= ##
## 21SPONGE Claire 30 Sources
sp30sc  = md.read_info_sponge_30src(fname = '../../oh/sponge/sub_data/30src_claire.txt')
sc1     = sp30sc['src']
xl1     = sp30sc['l']
xb1     = sp30sc['b']
hi1     = sp30sc['nhi']
hi1er   = sp30sc['nhi_er']
thin1   = sp30sc['thin']
thin1er = sp30sc['thin_er']
cnm1    = sp30sc['cnm']
cnm1er  = sp30sc['cnm_er']
wnm1    = sp30sc['wnm']
wnm1er  = sp30sc['wnm_er']
spdat   = read2dict(sp30sc)

## 78 MS Sources
# ms79sc  = md.read_info_ms_79sc(fname = '../result/nhi_lb_79src_HT03.txt') ## from HT03
ms79sc  = md.read_info_ms_79sc(fname = '../scale_NHI_HT03/nhi_79src_recal.txt') ## Recal. from HT03
sc2     = ms79sc['src']
xl2     = ms79sc['l']
xb2     = ms79sc['b']
hi2     = ms79sc['nhi']
hi2er   = ms79sc['nhi_er']
thin2   = ms79sc['thin']
thin2er = ms79sc['thin_er']
cnm2    = ms79sc['cnm']
cnm2er  = ms79sc['cnm_er']
wnm2    = ms79sc['wnm']
wnm2er  = ms79sc['wnm_er']
msdat   = read2dict(ms79sc)

print 'MS length:'
print len(sc2)

print 'SP length:'
print len(sc1)

print 'Merged-list length:'
mergedlist = list(set(sc1 + sc2))
print len(mergedlist)

## Common & different Sources
comsc = []  ## 14 src
difsc = []  ## 16 src
for sc in spdat:
	if (sc in msdat):
		comsc.append(sc)
	else:
		difsc.append(sc)

k = 0
for i in range(len(sc2)):
	k  = k + 1 
	sc = sc2[i]
	if(sc not in comsc):
		print ('{}    {}\t{:08.4f}  {:08.4f}  {:06.2f}  {:06.2f}  {:08.4f}  {:08.4f}  {:06.2f}  {:08.4f}  {:06.2f}  {:08.4f}'\
			.format(k, sc, xl2[i], xb2[i], hi2[i], hi2er[i], thin2[i], thin2er[i], cnm2[i], cnm2er[i], wnm2[i], wnm2er[i]  ))         ## nhi_thin_cnm_wnm_93src.txt
	else:
		j = sc1.index(sc)
		print ('{}    {}\t{:08.4f}  {:08.4f}  {:06.2f}  {:06.2f}  {:08.4f}  {:08.4f}  {:06.2f}  {:08.4f}  {:06.2f}  {:08.4f}'\
			.format(k, sc, xl1[j], xb1[j], hi1[j], hi1er[j], thin1[j], thin1er[j], cnm1[j], cnm1er[j], wnm1[j], wnm1er[j]  ))         ## nhi_thin_cnm_wnm_93src_sponge_prior.txt

for i in range(len(difsc)):
	k  = k + 1
	sc = difsc[i]
	print ('{}    {}\t{:08.4f}  {:08.4f}  {:06.2f}  {:06.2f}  {:08.4f}  {:08.4f}  {:06.2f}  {:08.4f}  {:06.2f}  {:08.4f}'\
		.format(k, sc, spdat[sc]['l'], spdat[sc]['b'], spdat[sc]['nhi'], spdat[sc]['nhi_er'],\
		  spdat[sc]['thin'], spdat[sc]['thin_er'], spdat[sc]['cnm'], spdat[sc]['cnm_er'], spdat[sc]['wnm'], spdat[sc]['wnm_er'] )) ## nhi_thin_cnm_wnm_93src.txt