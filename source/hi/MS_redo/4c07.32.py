import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import module            as md

from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from gfitflex       import gfit

## Get the index of a given velocity #
 #
 # params list v-axis List of Velocity_axis
 # params float vel Velocity
 #
 # return int idx Index of vel in List of velocity_axis
 # 
 # Author Van Hiep ##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

## Get Vrange Indexes ##
 #
 # params 1-D-array v     VLSR
 # params float     lowv  Lower limit
 # params float     upv   Upper limit
 #
 # return list
 #
 # version 01/2017 
 # author Nguyen Van Hiep ##
def get_vrange_id(v, lowv, upv):
	vmax = get_vel_index(v, lowv)
	vmin = get_vel_index(v, upv)
	return [vmin, vmax]

## Multiple (N) Gaussians + offset. ##
 #
 # params list  v    VLSR
 # params float zr   estimated constant zero offset of the data points.
 # params list  h    the array of N estimated heights of the Gaussians.
 # params list  v0   the array of N estimated centers of the Gaussians.
 # params list  w    the array of N estimated halfwidths of the Gaussians.
 #
 # return 1-D-array  tf  The calculated points.
 #
 # version 01/2017 
 # author Nguyen Van Hiep ##
def gcurv(v, zr, h, v0, w):
	#DETERMINE NR OF GAUSSIANS...
	ng = len(h)

	v  = np.array(v)
	h  = np.array(h)
	v0 = np.array(v0)
	w  = np.array(w)

	tf = 0.*v + zr
	for i in range(ng):
		if (w[i] > 0.):
			tf = tf + h[i]*np.exp(- ( (v-v0[i])/(0.6005612*w[i]))**2) # 0.6005612 - 1/e width

	return tf


## ============= MAIN ================ ##
## Class
fit  = gfit()

src  = '4C07.32'
print 'Fitting...' + src

specs = md.read_hi_specs(fname = '../data/nhi_opac_specs.txt')
v     = specs[src]['v']
v     = np.array(v)

Texp  = specs[src]['Texp']
Texp  = np.array(Texp)

tau   = specs[src]['tau']
tau   = np.array(tau)

## ABS
vlsr    = v
tau     = tau
# sigtau  = data.sigma
cont    = 0.

## EM - Set the velocity range for the emission spectra (to trim Carl's data)
vlsrem  = v
vmin, \
vmax    = get_vrange_id(vlsrem, -100., 60.)
Te      = Texp
xdataem = vlsrem[vmin:vmax+1]
Te      = Te[vmin:vmax+1]

## Absorption line ###
plt.figure(1, figsize=(12,12))
plt.plot(vlsr,tau, 'k-', linewidth=2, label='data, Absorption line')
plt.ylabel(r'$\tau$', fontsize=35)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper left', fontsize=18)
plt.grid(False)
plt.xlim(-50, 50)
plt.show()

## Emission lines ###
plt.figure(2, figsize=(12,12))
plt.plot(xdataem, Te, 'k-', linewidth=2, label='data, Emission line')
plt.ylabel(r'$Texp$', fontsize=35)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(-50, 50)
plt.show()

# Retrieve initial Gaussian guesses for the absorption spectrum
zro0   = 0.0
hgt0   = [0.177]
cen0   = [-1.9]
wid0   = [2.19]
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

plt.figure(3, figsize=(12,12))
plt.plot(vlsr, tfita, 'r-', linewidth=2, label='fit')
plt.plot(vlsr,tau, 'b-', linewidth=2, label='data, Absorption line')
plt.ylabel(r'$\tau$', fontsize=35)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper left', fontsize=18)
plt.grid(False)
plt.xlim(-50, 50)
plt.show()



zrownm1   = 0.0
hgtwnm1   = [0.005]
cenwnm1   = [-10.2]
widwnm1   = [30.49]
look      = 0
nrgwnm    = len(hgtwnm1)
zrownm1yn = 1
hgtwnm1yn = [1]*nrgwnm
cenwnm1yn = [1]*nrgwnm
widwnm1yn = [1]*nrgwnm
fwnm1     = [0]*nrgwnm
fwnm1yn   = [0]*nrgwnm

order1    = list(range(nrg))
tspin1    = [50.]*nrg
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
tsky  = 3.3711
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

print zrocnm1, hgtcnm1, cencnm1, widcnm1, tspincnm1
print 'Expected Emission profile: problem...', problem

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


plt.figure(4, figsize=(12,12))
plt.plot(xdataem, tdata, 'k-', linewidth=2, label='data, Emission line')
plt.plot(xdataem, tb_tot_fit, 'r-', label='Fit, Emission line')
plt.ylabel(r'$Texp$', fontsize=35)
plt.tick_params(axis='y', labelsize=15)
plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(-50, 50)
plt.show()


sys.exit()