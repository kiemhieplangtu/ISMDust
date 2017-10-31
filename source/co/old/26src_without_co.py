import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import pylab             as pl
import module            as md

from numpy    import array
from restore  import restore

## ========= MAIN ============= ##
inf26 = md.read_basic_info_no_co(fname = 'result/26src_no_co_basic.dat')
sc26  = inf26['src']
l     = inf26['l']
b     = inf26['b']
rai   = inf26['ra_icrs']
deci  = inf26['de_icrs']
raj   = inf26['ra_j']
decj  = inf26['de_j']
oh    = inf26['oh']

## 94 HI sources
inf94  = md.read_nhi_93src(fname = '../hi/result/nhi_thin_cnm_wnm_93src_scaled.txt')
sc94   = inf94['src']
nhi    = inf94['nhi']
nhier  = inf94['nhi_er']
thin   = inf94['thin']
thiner = inf94['thin_er']
cnm    = inf94['cnm']
cnmer  = inf94['cnm_er']
wnm    = inf94['wnm']
wnmer  = inf94['wnm_er']

for i in range(len(sc26)):
	k = sc94.index(sc26[i])
	print('{:3}  {:11} {:.4f}    {:.4f}    {:.4f}    {:.4f}    {}    {}    {}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}    {:.4f}'\
		.format(i, sc26[i], l[i], b[i], rai[i], deci[i], raj[i], decj[i], oh[i], nhi[k], nhier[k], thin[k], thiner[k], cnm[k], cnmer[k], wnm[k], wnmer[k] ) )      ## result/26src_no_co_with_sponge.dat