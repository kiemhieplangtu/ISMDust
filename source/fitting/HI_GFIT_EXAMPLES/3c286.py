import sys, os
sys.path.insert(0, 'lib') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import module            as md

from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from gfitflex       import gfit

## Read info of HI EM and ABS spectra #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 09/2017
 # Author Van Hiep ##
def read_hi_specs_3c286(fname = 'data/3c286_ABS_data.txt'):
	cols  = ['id', 'src','v', 'yy']
	fmt   = ['i',  's',  'f', 'f'  ]
	data  = restore(fname, 2, cols, fmt)
	dat   = data.read(asarray=True)

	src   = dat['src']
	v     = dat['v']
	yy    = dat['yy']

	return v,yy


## ============= MAIN ================ ##
## Class
fit  = gfit()

src  = '3C286'
print 'Fitting...' + src

# data = readsav('gfit_idl_claire/3C286.sav')
# % RESTORE: Restored variable: VLSR.
# % RESTORE: Restored variable: SPEC1.
# % RESTORE: Restored variable: EM_SPEC.
# % RESTORE: Restored variable: SIGMA.
# % RESTORE: Restored variable: EMSIGMA.
# % RESTORE: Restored variable: PSR1.
# % RESTORE: Restored variable: VLSREM.
# % RESTORE: Restored variable: CONT.
# % RESTORE: Restored variable: RMS.

vlsr,tau  = read_hi_specs_3c286(fname = 'data/3c286_ABS_data.txt')
vlsrem,Te = read_hi_specs_3c286(fname = 'data/3c286_EM_data.txt')

## ABS
# vlsr    = data.vlsr[:,0]
# tau     = data.spec1[:,0,0]
# sigtau  = data.sigma
# cont    = data.cont
cont    = 0.

# print vlsr
# print tau
# print sigtau
# print cont

## EM - Set the velocity range for the emission spectra (to trim Carl's data)
# vlsrem  = data.vlsrem
vmin, \
vmax    = md.get_vrange_id(vlsrem, -100., 100.)
# Te      = data.em_spec	

xdataem = vlsrem[vmin:vmax+1]
Te      = Te[vmin:vmax+1]


## Absorption profile ###
plt.figure(1, figsize=(12,12))
plt.plot(vlsr,tau, 'k-', linewidth=2, label='data, Absorption line')
plt.title(src)
plt.xlabel('Vlsr', fontsize=35)
plt.ylabel(r'$\tau$', fontsize=35)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper left', fontsize=18)
plt.grid(False)
plt.xlim(-50, 50)
plt.show()

## Emission profile ###
plt.figure(2, figsize=(12,12))
plt.plot(xdataem, Te, 'k-', linewidth=2, label='data, Emission line')
plt.title(src)
plt.xlabel('Vlsr', fontsize=35)
plt.ylabel('Texp', fontsize=35)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(-50, 50)
plt.show()



### Do CNM ###
# Retrieve initial Gaussian guesses for the absorption spectrum
zro0   = 0.0
hgt0   = [0.0065, 0.0066, 0.009]
cen0   = [-28.48, -14.4, -7.3]
wid0   = [2.32,    4.66,   5.0]

hgt0   = [1., 1., 1.]
cen0   = [-30., -15., -5.]
wid0   = [1., 1., 1.]
look   = 0

nrg    = len(hgt0)
zro0yn = 0
hgt0yn = [1]*nrg
cen0yn = [1]*nrg
wid0yn = [1]*nrg
corder = 'no'

## WNM
tspin0 = [0.]*nrg
order0 = list(range(nrg))

zro0yn   = 0
tau0yn   = [1]*nrg
cenc0yn  = [1]*nrg
wid0yn   = [1]*nrg
tspin0yn = [0]*nrg

zrownm = 1.
hgtwnm = [0.003]
cenwnm = [-10.4]
widwnm = [21.46]
fwnm   = [0.5]

zrownmyn = 1
hgtwnmyn = 0
cenwnmyn = 0
widwnmyn = 0
fwnmyn   = 0

## Fit these guesses
tfita, sigma, resida,\
zro1, hgt1, cen1, wid1,\
sigzro1, sighgt1, sigcen1, sigwid1,\
cov, problem,\
nparams = fit.fit(look, vlsr, tau, [0, len(tau)-1],\
			      zro0, hgt0, cen0, wid0,\
				  zro0yn, hgt0yn, cen0yn, wid0yn)

print 'Absorption line: problem...', problem

print '1. sigma ', sigma
print '2. Zro ', zro1
print '3. tau ', hgt1
print '\t\t', sighgt1
print '4. cen ', cen1
print '\t\t', sigcen1
print '5. wid ', wid1
print '\t\t', sigwid1

print ''
print ''

### Plot: Absorption data and Fit ###
plt.figure(3, figsize=(12,12))
plt.plot(vlsr, tfita, 'r-', linewidth=2, label='fit')
plt.plot(vlsr,tau, 'k-', linewidth=1, label='data, Absorption line')
plt.title(src)
plt.xlabel('Vlsr', fontsize=35)
plt.ylabel(r'$\tau$', fontsize=35)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper left', fontsize=18)
plt.grid(False)
plt.xlim(-50, 50)
plt.show()







### Do WNM ##
zrownm1   = 0.0
hgtwnm1   = [1,1]
cenwnm1   = [-5,-20]
widwnm1   = [10,10]
look      = 0
nrgwnm    = len(hgtwnm1)
zrownm1yn = 1
hgtwnm1yn = [1]*nrgwnm
cenwnm1yn = [1]*nrgwnm
widwnm1yn = [1]*nrgwnm
fwnm1     = [0]*nrgwnm
fwnm1yn   = [0]*nrgwnm

order1    = list(range(nrg))
tspin1    = [30.]*nrg
tspin1yn  = [1]*nrg

zrocnm1   = 0.
hgtcnm1   = hgt1
cencnm1   = cen1
widcnm1   = wid1
zrocnm1yn = 0
hgtcnm1yn = [0]*nrg
cencnm1yn = [0]*nrg
widcnm1yn = [0]*nrg

look      = -1
xindxrange= [0,len(xdataem)-1]

## ---Parameters within tbgfitflex_exp.pro, sets number of loops (nloopmax)
## ---and the fractional change in each parameter per loop iteration
nloopmax     = 100
halfasseduse = 0.2

## ---Compute Tsky at the source position [**predict_sync method currently not working, 
## ---so I just set a generic value here, will fix in the future]
## @galactic_coordinates.pro
## print, 'gl gb', gl, ' ', gb
## tsky=predict_sync(gl,gb, nu=1.4, /interp)+2.725 
tsky  = 2.8 ## 3.41
tdata = Te+tsky

tfite, sigma, reside,\
zrocnm1, hgtcnm1, cencnm1, widcnm1, tspincnm1, \
sigzrocnm1, sighgtcnm1, sigcencnm1, sigwidcnm1, sigtspincnm1, \
zrownm1, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
sigzrownm1, sighgtwnm1, sigcenwnm1, sigwidwnm1, sigfwnm1, \
cov, problem, nloop, \
tb_cont, tb_wnm_tot, tb_cnm_tot, \
exp_tau_sum, nloopmax, halfasseduse = fit.efit(look, xdataem, tdata, xindxrange, \
			zrocnm1, hgtcnm1, cencnm1, widcnm1, tspin1, order1, \
			zrocnm1yn, hgtcnm1yn, cencnm1yn, widcnm1yn, tspin1yn, \
			cont, hgtwnm1, cenwnm1, widwnm1, fwnm1, \
			1, hgtwnm1yn, cenwnm1yn, widwnm1yn, fwnm1yn, nloopmax=nloopmax, halfasseduse=0.2)

tb_tot_fit = tb_cnm_tot+tb_wnm_tot+tb_cont-tsky
tb_tot_fit = tfite-tsky
tdata      = tdata-tsky

# print zrocnm1, hgtcnm1, cencnm1, widcnm1, tspincnm1
print 'Expected Emission profile: problem...', problem
print 'Nloop:', nloop
print 'Nloop-max', nloopmax

print '1. sigma ', sigma
print '2. Zro ', zro1
print '3. tau ', hgt1
print '\t\t', sighgt1
print '4. cen ', cen1
print '\t\t', sigcen1
print '5. wid ', wid1
print '\t\t', sigwid1
print '6. Tspin ', tspincnm1
print '\t\t', sigtspincnm1

print ''
print ''

print '7. Tau-WNM ', hgtwnm1
print '\t\t', sighgtwnm1
print '8. V0-WNM ', cenwnm1
print '\t\t', sigcenwnm1
print '9. Width-WNM ', widwnm1
print '\t\t', sigwidwnm1
print '10. fwnm1-WNM ', fwnm1
print '\t\t', sigfwnm1


### Plot: Emission data and Fit ###
plt.figure(4, figsize=(12,12))
plt.plot(xdataem, tdata, 'k-', linewidth=2, label='data, Emission line')
plt.plot(xdataem, tb_tot_fit, 'r-', label='Fit, Emission line')
plt.title(src)
plt.xlabel('Vlsr', fontsize=35)
plt.ylabel('Texp', fontsize=35)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(-50, 50)
plt.show()


sys.exit()
### END ###