import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import healpy            as hp
import module            as md
import copy

from restore             import restore


## Read info of OH sources #
 # l,b, noh, noh_er
 #
 # params string fname Filename
 # return dict info
 # 
 # version 4/2017
 # Author Van Hiep ##
def read_total_noh(fname = '../../oh/result/total_noh65_21src.txt'):
	cols = ['idx','src','l', 'b', 'noh', 'noh_er']
	fmt  = ['i', 's',   'f', 'f', 'f',   'f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	noh  = dat['noh']
	er   = dat['noh_er']
	src  = dat['src']
	l    = dat['l']
	b    = dat['b']

	ret  = {}
	for i in range(0,len(src)):
		# if dat['src'][i] not in ret.keys():
		ret[src[i]] = {}
		ret[src[i]]['noh']   = noh[i]
		ret[src[i]]['l']     = l[i]
		ret[src[i]]['b']     = b[i]
		ret[src[i]]['noher'] = er[i]

	return ret

##================= MAIN ========================##
info   = md.read_nhi_93src(fname = '../../hi/result/nhi_thin_cnm_wnm_94src_scaled.txt')
noh    = read_total_noh(fname = 'total_noh67_21src.txt')
# noh  = read_noh(fname = '../../oh/NOH/total_noh67_21src_carl.txt')

src    = info['src']
xl     = info['l']
xb     = info['b']
nhi    = info['nhi']
nhier  = info['nhi_er']
thin   = info['thin']
thiner = info['thin_er']
cnm    = info['cnm']
cnmer  = info['cnm_er']
wnm    = info['wnm']
wnmer  = info['wnm_er'] 

print src, len(src)

k = 0
for i in range(len(src)):
	if(src[i] not in noh):
		continue

	print k, src[i], xl[i], xb[i], nhi[i], nhier[i], thin[i], thiner[i], cnm[i], cnmer[i], wnm[i], wnmer[i], noh[ src[i] ]['l'], noh[ src[i] ]['b'], noh[ src[i] ]['noh'], noh[ src[i] ]['noher']
	k = k + 1 