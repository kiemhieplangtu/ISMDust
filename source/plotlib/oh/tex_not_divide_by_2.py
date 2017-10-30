import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import matplotlib
import pandas            as pd

from scipy.optimize  import curve_fit
from scipy.io.idl    import readsav
from numpy           import array
from restore         import restore
from plotting        import cplot
from gauss_fit       import gfit
from scipy.integrate import quad

## Create tau function ##
 #
 # params list x x axis
 # params list params Parameters
 #
 # paurn array y Multiple-Gaussian functions
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def tau_func(x, *params):
	y = 0.
	for i in range(0, len(params), 4):
		amp = params[i]
		ctr = params[i+1]
		wid = params[i+2]
		wid = wid/2./np.sqrt(2.*np.log(2.))
		tex = params[i+3]
		# y   = y + tex*amp * np.exp( -0.5*((x - ctr)/wid)**2 )
		y   = y + amp * np.exp( -0.5*((x - ctr)/wid)**2 )
	return y

#============== MAIN ==============#
dat                 = readsav('data/doall_fitgauss.sav')

## For OH1665 Absorption 
src1                = dat.amg_1665.sname
ell1                = dat.amg_1665.ell
bee1                = dat.amg_1665.bee
ng1                 = dat.amg_1665.nrgauss
gaussnr1            = dat.amg_1665.gaussnr
zro1                = dat.amg_1665.zrocnm
taucnm1             = dat.amg_1665.taucnm
cencnm1             = dat.amg_1665.cencnm
widcnm1             = dat.amg_1665.widcnm
tspincnm1           = dat.amg_1665.tspincnm
sigzrocnm1          = dat.amg_1665.sigzrocnm
sigtaucnm1          = dat.amg_1665.sigtaucnm
sigcencnm1          = dat.amg_1665.sigcencnm
sigwidcnm1          = dat.amg_1665.sigwidcnm
sigtspincnm1        = dat.amg_1665.sigtspincnm
continuum_defln1    = dat.amg_1665.continuum_defln
sigcontinuum_defln1 = dat.amg_1665.sigcontinuum_defln
continuum_em1       = dat.amg_1665.continuum_em


## For OH1667 Absorption 
src2                = dat.amg_1667.sname
ng2                 = dat.amg_1667.nrgauss
gaussnr2            = dat.amg_1667.gaussnr
ell2                = dat.amg_1667.ell
bee2                = dat.amg_1667.bee
zro2                = dat.amg_1667.zrocnm
taucnm2             = dat.amg_1667.taucnm
cencnm2             = dat.amg_1667.cencnm
widcnm2             = dat.amg_1667.widcnm
tspincnm2           = dat.amg_1667.tspincnm
sigzrocnm2          = dat.amg_1667.sigzrocnm
sigtaucnm2          = dat.amg_1667.sigtaucnm
sigcencnm2          = dat.amg_1667.sigcencnm
sigwidcnm2          = dat.amg_1667.sigwidcnm
sigtspincnm2        = dat.amg_1667.sigtspincnm
continuum_defln2    = dat.amg_1667.continuum_defln
sigcontinuum_defln2 = dat.amg_1667.sigcontinuum_defln
continuum_em2       = dat.amg_1667.continuum_em

## For OH1665 Emmission
esrc1                = dat.emg_1665.sname
# elle1                = dat.emg_1665.ell
# beee1                = dat.emg_1665.bee
eng1                 = dat.emg_1665.nrgauss
bsl1                 = dat.emg_1665.tbaseline
cont_em1             = dat.emg_1665.continuum_em
hgtwnm1              = dat.emg_1665.hgtwnm
cenwnm1              = dat.emg_1665.cenwnm
widwnm1              = dat.emg_1665.widwnm
fwnm1                = dat.emg_1665.fwnm
sigcont1             = dat.emg_1665.sigcontinuum
sighgtwnm1           = dat.emg_1665.sighgtwnm
sigcenwnm1           = dat.emg_1665.sigcenwnm
sigwidwnm1           = dat.emg_1665.sigwidwnm
sigfwnm1             = dat.emg_1665.sigfwnm

## For OH1667 Emmission 
esrc2                = dat.emg_1667.sname
# elle2                = dat.emg_1667.ell
# beee2                = dat.emg_1667.bee
eng2                 = dat.emg_1667.nrgauss
bsl2                 = dat.emg_1667.tbaseline
cont_em2             = dat.emg_1667.continuum_em
hgtwnm2              = dat.emg_1667.hgtwnm
cenwnm2              = dat.emg_1667.cenwnm
widwnm2              = dat.emg_1667.widwnm
fwnm2                = dat.emg_1667.fwnm
sigcont2             = dat.emg_1667.sigcontinuum
sighgtwnm2           = dat.emg_1667.sighgtwnm
sigcenwnm2           = dat.emg_1667.sigcenwnm
sigwidwnm2           = dat.emg_1667.sigwidwnm
sigfwnm2             = dat.emg_1667.sigfwnm

print esrc2
print esrc1
print hgtwnm2
for i in range(len(esrc2)):
	if(esrc2[i] in esrc1):
		if(i>0):
			j=i-1
		else:
			j=i
		print esrc2[i], '&', str(round(ell2[i],1) )+ '/'+ str(round(bee2[i],1) ), '&', \
		round(hgtwnm1[j],4),'$\pm$', round(sighgtwnm1[j],4), '&', \
		round(cenwnm1[j],2),'$\pm$', round(sigcenwnm1[j],2), '&', \
		round(widwnm1[j],2),'$\pm$',round(sigwidwnm1[j],2), '&', \
		'- &', \
		'- &', round(fwnm1[j],2),'$\pm$',round(sigfwnm1[j],2), '& & ', \
		round(hgtwnm2[i],4),'$\pm$', round(sighgtwnm2[i],4), '&', \
		round(cenwnm2[i],2),'$\pm$', round(sigcenwnm2[i],2), '&', \
		round(widwnm2[i],2),'$\pm$',round(sigwidwnm2[i],2), '&', \
		'- &', \
		'- &', round(fwnm2[i],2),'$\pm$',round(sigfwnm2[i],2), '\\\\'  ## oh_infor_to_latex_table.txt
	else:
		print esrc2[i], '&', str(round(ell2[i],1) )+ '/'+ str(round(bee2[i],1) ), '&', \
		'- &', \
		'- &', \
		'- &', \
		'- &', \
		'- &', '-', '& & ', \
		round(hgtwnm2[i],4),'$\pm$', round(sighgtwnm2[i],4), '&', \
		round(cenwnm2[i],2),'$\pm$', round(sigcenwnm2[i],2), '&', \
		round(widwnm2[i],2),'$\pm$',round(sigwidwnm2[i],2), '&', \
		'- &', \
		'-&', round(fwnm2[i],2),'$\pm$',round(sigfwnm2[i],2), '\\\\'  ## oh_infor_to_latex_table.txt

c2  = 2.2230287734  #e14
c1  = 3.99757843817 #e14
pa1 = {}
pa2 = {}
p1  = {}
p2  = {}
for i in range(len(taucnm1)):
	p1[i]  = [taucnm1[i], cencnm1[i], widcnm1[i],tspincnm1[i]]
	p2[i]  = [taucnm2[i], cencnm2[i], widcnm2[i],tspincnm2[i]]
	if src2[i] not in pa2.keys():
		k               = 0
		pa2[src2[i]]    = {}
		pa2[src2[i]][k] =  [taucnm2[i], cencnm2[i], widcnm2[i],tspincnm2[i]]
	else:
		k               = k + 1
		pa2[src2[i]][k] = [taucnm2[i], cencnm2[i], widcnm2[i],tspincnm2[i]]

	if src1[i] not in pa1.keys():
		j               = 0
		pa1[src1[i]]    = {}
		pa1[src1[i]][j] =  [taucnm1[i], cencnm1[i], widcnm1[i],tspincnm1[i]]
	else:
		j               = j + 1
		pa1[src1[i]][j] = [taucnm1[i], cencnm1[i], widcnm1[i],tspincnm1[i]]


## N(OH)
pi   = np.pi
xmin = -100.0
xmax = +100.0
for i in range(len(taucnm1)):
	par2  = p2[i]
	stau  = quad(tau_func, xmin,xmax, args=tuple(par2))
	noh67 = c2 * par2[3] * stau[0]  # x10^14 
	noh2  = 2.36634 * par2[0]* par2[3] * par2[2] # x10^14
	sta   = par2[0]*par2[2]*np.sqrt(pi)/2./np.sqrt(np.log(2.))
	if(par2[3] == 1.0):
		noh67 = 0.
		noh2  = 0.
	
	par1  = p1[i]
	stau  = quad(tau_func, xmin,xmax, args=tuple(par1))
	noh65 = c1 * par1[3] * stau[0]  # x10^14 
	noh1  = 4.26413 * par1[0]* par1[3] * par1[2] # x10^14
	sta   = par1[0]*par1[2]*np.sqrt(pi)/2./np.sqrt(np.log(2.))
	if(par1[3] == 1.0):
		noh65 = 0.
		noh1  = 0.
	# print src2[i], '   \t', str(noh2)+'/'+ str(noh67), '   \t', str(noh1)+'/'+ str(noh65)

## N(OH) - Simpler
noh2 = 2.36634* taucnm2 * tspincnm2 * widcnm2
noh1 = 4.26413* taucnm1 * tspincnm1 * widcnm1

dd1  = widcnm1**2   * sigtspincnm1**2 * taucnm1 + \
       widcnm1**2   * tspincnm1**2    * sigtaucnm1**2 + \
       tspincnm1**2 * taucnm1**2      * sigwidcnm1**2
d1   = 4.26413*np.sqrt(dd1)

dd2  = widcnm2**2   * sigtspincnm2**2 * taucnm2 + \
       widcnm2**2   * tspincnm2**2    * sigtaucnm2**2 + \
       tspincnm2**2 * taucnm2**2      * sigwidcnm2**2
d2   = 2.36634*np.sqrt(dd2)

# for i in range(len(taucnm1)):
# 	if(tspincnm2[i] == 1.0):
# 		noh2[i]  = 0.
# 	if(tspincnm1[i] == 1.0):
# 		noh1[i]  = 0.
# 	print src2[i], '   \t', str(noh2[i])+'/'+ str(d2[i]), '   \t', str(noh1[i])+'/'+ str(d1[i])


# for i in range(len(src1)):
# 	print src1[i], '&', str(round(ell1[i],1) )+ '/'+ str(round(bee1[i],1) ), '&', \
# 	round(taucnm1[i],4),'$\pm$', round(sigtaucnm1[i],4), '&', \
# 	round(cencnm1[i],2),'$\pm$', round(sigcencnm1[i],2), '&', \
# 	round(widcnm1[i],2),'$\pm$',round(sigwidcnm1[i],2), '&', \
# 	round(tspincnm1[i],2),'$\pm$',round(sigtspincnm1[i],2), '&', \
# 	round(noh1[i],2),'$\pm$', round(d1[i],2), '&', gaussnr1[i], '& & ', \
# 	round(taucnm2[i],4),'$\pm$', round(sigtaucnm2[i],4), '&', \
# 	round(cencnm2[i],2),'$\pm$', round(sigcencnm2[i],2), '&', \
# 	round(widcnm2[i],2),'$\pm$',round(sigwidcnm2[i],2), '&', \
# 	round(tspincnm2[i],2),'$\pm$',round(sigtspincnm2[i],2), '&', \
# 	round(noh2[i],2),'$\pm$', round(d2[i],2), '&', gaussnr2[i], '\\\\'  ## oh_infor_to_latex_table.txt

# for i in range(len(src1)):
# 	print i, src1[i], '\t', round(ell1[i],3),'\t', round(bee1[i],3), '\t', \
# 	round(taucnm1[i],7),'\t', round(sigtaucnm1[i],7), '\t', \
# 	round(cencnm1[i],3),'\t', round(sigcencnm1[i],3), '\t', \
# 	round(widcnm1[i],3),'\t',round(sigwidcnm1[i],3), '\t', \
# 	round(tspincnm1[i],3),'\t',round(sigtspincnm1[i],3), '\t', \
# 	round(noh1[i],4),'\t', round(d1[i],4), '\t', gaussnr1[i], '  ', \
# 	round(taucnm2[i],7),'\t', round(sigtaucnm2[i],7), '\t', \
# 	round(cencnm2[i],3),'\t', round(sigcencnm2[i],3), '\t', \
# 	round(widcnm2[i],3),'\t',round(sigwidcnm2[i],3), '\t', \
# 	round(tspincnm2[i],3),'\t',round(sigtspincnm2[i],3), '\t', \
# 	round(noh2[i],4),'\t', round(d2[i],4), '\t', gaussnr2[i]

for i in range(len(src1)):
	print i, src1[i], '\t', str(round(ell1[i],1) )+ '/'+ str(round(bee1[i],1) ), '\t', \
	round(taucnm1[i],4),'/', round(sigtaucnm1[i],4), '\t', \
	round(cencnm1[i],2),'/', round(sigcencnm1[i],2), '\t', \
	round(widcnm1[i],2),'/',round(sigwidcnm1[i],2), '\t', \
	round(tspincnm1[i],2),'/',round(sigtspincnm1[i],2), '\t', \
	round(noh1[i],2),'/', round(d1[i],2), '  ', \
	round(taucnm2[i],4),'/', round(sigtaucnm2[i],4), '\t', \
	round(cencnm2[i],2),'/', round(sigcencnm2[i],2), '\t', \
	round(widcnm2[i],2),'/',round(sigwidcnm2[i],2), '\t', \
	round(tspincnm2[i],2),'/',round(sigtspincnm2[i],2), '\t', \
	round(noh2[i],2),'/', round(d2[i],2)