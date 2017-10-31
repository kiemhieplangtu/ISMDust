import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import operator

from numpy    import array
from restore  import restore
from plotting import cplot

## Compute N(OH) ##
 #
 # params dict dat Data on N(OH)
 # return dict ret N(OH) of each source
 #
 # version 09/2016 
 # author Nguyen Van Hiep ##
def get_noh_for_each_src(fname, skip=3):
	cols = ['src', 'l', 'b', 'tau1', 'tau1er', 'v01', 'v01er','wid1', 'wid1er','tex1','tex1_er','noh1','noh1er', 'tau2', 'tau2er', 'v02', 'v02er','wid2', 'wid2er','tex2','tex2_er','noh2','noh2er']
	fmt  = ['s',   'f', 'f', 'f',    'f',      'f',    'f',   'f',     'f',     'f',  'f',       'f',   'f',       'f',    'f',      'f',    'f',   'f',     'f',     'f',  'f',       'f',   'f'    ]
	data = restore(fname, skip, cols, fmt)
	dat  = data.read()
	noh  = dat['noh2']
	er2  = dat['noh2er']
	src  = dat['src']
	l    = dat['l']
	b    = dat['b']

	ret  = {}
	for i in range(0,len(src)):
		if dat['src'][i] not in ret.keys():
			ret[src[i]]        = {}
			ret[src[i]]['noh'] = noh[i]
			ret[src[i]]['l']   = l[i]
			ret[src[i]]['b']   = b[i]
			ret[src[i]]['er2'] = er2[i]**2
		else:
			ret[src[i]]['l']   = l[i]
			ret[src[i]]['b']   = b[i]
			ret[src[i]]['noh'] = ret[src[i]]['noh'] + noh[i]
			ret[src[i]]['er2'] = ret[src[i]]['er2'] + er2[i]**2

	return ret

##================= MAIN ========================##
## N(OH) ##
# noh_toff_binup.py > infor_oh_components_simple.txt
# noh   = get_noh_for_each_src(fname = 'infor_oh_components_simple.txt')
noh   = get_noh_for_each_src(fname = 'infor_oh_components_carl.txt', skip=4)
# noh   = get_noh_for_each_src(fname = 'noh_components_simple_raw.txt')
n     = 0

file  = open('total_raw_noh67_21src.txt','w') 
for sc in noh:
	# print n, '\t', sc, '\t', round(noh[sc]['l'],8), '\t', round(noh[sc]['b'],8), '\t', round(noh[sc]['noh'],8), '\t', round(np.sqrt(noh[sc]['er2']),8)
	text = str(n)+ '\t'+ sc+ '\t'+ str( round(noh[sc]['l'],8) )+ '\t' + str( round(noh[sc]['b'],8) )+ '\t' + str( round(noh[sc]['noh'],8) )+ '\t'+ str( round(np.sqrt(noh[sc]['er2']),8) ) + '\n'
	file.write( text )
	n    = n + 1

file.close()

## total_noh_21src.txt
## total_raw_noh_21src.txt