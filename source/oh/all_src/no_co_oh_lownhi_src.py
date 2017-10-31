import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import module            as md

from scipy.io.idl        import readsav
from restore             import restore
from mpfit               import mpfit


#============== MAIN ==============#
# Info of 26 sources with no CO - l/b/name && 23 src low NHI #
noco   = md.read_info_no_co('../../co12/result/26src_no_CO_scaled.txt')
lownhi = md.read_lownhi_23src(fname = '../../hi/result/lownhi_thin_cnm_wnm_scaled.txt')

## Print both of them, or one by one ##
## NO CO
sc     = noco['src']
xl     = noco['l']
xb     = noco['b']
hi     = noco['nhi']
hier   = noco['nhi_er']
thin   = noco['thin']
thiner = noco['thin_er']
cnm    = noco['cnm']
cnmer  = noco['cnm_er']
wnm    = noco['wnm']
wnmer  = noco['wnm_er']
ohyn   = noco['oh']

k = 0
for i in range(len(noco['nhi'])):
	if(ohyn[i] == 0):
		print ('{}    {}\t{:08.4f}  {:08.4f}  {:06.2f}  {:06.2f}  {:08.4f}  {:08.4f}  {:06.2f}  {:08.4f}  {:06.2f}  {:08.4f}'\
			.format(k, sc[i], xl[i], xb[i], hi[i], hier[i], thin[i], thiner[i], cnm[i], cnmer[i], wnm[i], wnmer[i]  ))         ## 36src_noCO_noOH_lowNHI.txt
		k = k + 1
		## 19src_noCO_noOH.txt

## Low N(HI)
sc     = lownhi['src']
xl     = lownhi['l']
xb     = lownhi['b']
hi     = lownhi['nhi']
hier   = lownhi['nhi_er']
thin   = lownhi['thin']
thiner = lownhi['thin_er']
cnm    = lownhi['cnm']
cnmer  = lownhi['cnm_er']
wnm    = lownhi['wnm']
wnmer  = lownhi['wnm_er']
ohyn   = lownhi['oh']

k = 0
for i in range(len(lownhi['nhi'])):
	# if(ohyn[i] == 0):
	print ('{}    {}\t{:08.4f}  {:08.4f}  {:06.2f}  {:06.2f}  {:08.4f}  {:08.4f}  {:06.2f}  {:08.4f}  {:06.2f}  {:08.4f}'\
		.format(k, sc[i], xl[i], xb[i], hi[i], hier[i], thin[i], thiner[i], cnm[i], cnmer[i], wnm[i], wnmer[i]  ))         ## 36src_noCO_noOH_lowNHI.txt
	k = k + 1
	## 23src_lowNHI.txt

sys.exit()