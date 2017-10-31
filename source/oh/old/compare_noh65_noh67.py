import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import operator
import copy

from numpy               import array
from restore             import restore
from plotting            import cplot
from mpfit               import mpfit
from scipy.odr           import *

## Linear fucntion ##
 #
 # params list/float p Parameters
 # params 
 #
 # return 
 #
 # version 03/2017 
 # author Nguyen Van Hiep ##
def lin_fc(p, x):
     m = p
     return m*x

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

	return noh, er

## Read infor for each OH src
 # Fit inf from Carl's paper
 #
 # inf void
 #
 # return dict Infor for each src
 # 
 # version 2/2017
 # Author Van Hiep ##
def sum_noh(fname='oh_components_infor.txt'):
	## Read fit inf of HI components from Carl papers ##
	cols = ['id', 'src', 'l', 'b', 'tau1', 'tau1_er', 'v01', 'v01_er', 'wid1', 'wid1_er', 'tex1', 'tex1_er', 'noh1', 'noh1_er', 'ord1', 'tau2', 'tau2_er', 'v02', 'v02_er', 'wid2', 'wid2_er', 'tex2', 'tex2_er', 'noh2', 'noh2_er', 'ord2']
	fmt  = ['i',  's',  'f',  'f',  'f',   'f',        'f',  'f',      'f',    'f',       'f',     'f',      'f',     'f',      'i',     'f',    'f',       'f',   'f',      'f',    'f',       'f',     'f',      'f',    'f',       'i']
	dat  = restore(fname, 3, cols, fmt)
	inf  = dat.read()

	sc    = inf['src']
	noh1  = {}
	noh1e = {}

	noh2  = {}
	noh2e = {}
	for i in range(len(sc)):
		if sc[i] not in noh1.keys():
			noh1[sc[i]]  = {}
			noh1e[sc[i]] = {}

			noh2[sc[i]]  = {}
			noh2e[sc[i]] = {}

			noh1[sc[i]]  = inf['noh1'][i]
			noh1e[sc[i]] = inf['noh1_er'][i]

			noh2[sc[i]]  = inf['noh2'][i]
			noh2e[sc[i]] = inf['noh2_er'][i]
		else:
			noh1[sc[i]]  += inf['noh1'][i]
			noh1e[sc[i]] += inf['noh1_er'][i]

			noh2[sc[i]]  += inf['noh2'][i]
			noh2e[sc[i]] += inf['noh2_er'][i]

	oh1    = []
	oh2    = []
	oh1_er = []
	oh2_er = []
	for s in noh1:
		oh1.append(noh1[s])
		oh2.append(noh2[s])
		oh1_er.append(noh1e[s])
		oh2_er.append(noh2e[s])

	return oh1, oh2, oh1_er, oh2_er

##================= MAIN ========================##
## N(OH) from Heiles Method ##
xnoh1, xnoh2, xnoh1e, xnoh2e = sum_noh(fname='../plotlib/oh/oh_components_infor.txt')

## N(OH) from Simple method ##
noh1, noh1e = read_total_noh(fname = 'result/total_noh65_21src.txt')
noh2, noh2e = read_total_noh(fname = 'result/total_noh67_21src.txt')

## ********* Compare N(OH)65 and N(OH)67 from Simple method ********* ##
xdata = np.array(noh1)
ydata = np.array(noh2)

# Error bar for x-axis and y-axis
xerr  = np.array(noh1e)
yerr  = np.array(noh2e)

## Fit ##
lguess  = [ 1.]
########### ODR fit ############
# Create a model for Orthogonal distance regression (ODR) fitting.
lin_model = Model(lin_fc)
# Create a RealData object using our initiated data from above.
data      = RealData(xdata, ydata, sx=xerr, sy=yerr)
# Set up ODR with the model and data.
odr       = ODR(data, lin_model, beta0=lguess)
# Run the regression.
out       = odr.run()

print '********* ODR Results *********'
abp   = out.beta
abper = out.sd_beta
for i in range(1):
	print "%s = %03.8f +/- %03.8f" % ('slope',abp[i],abper[i])
## Plot ##
a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
# b     = np.array([ abp[1]-abper[1], abp[1]+abper[1] ])
xfit  = np.linspace(xdata.min(), xdata.max(), 20)
# yfit  = a[:, None] * xfit + b[:, None]
yfit  = a[:, None] * xfit
mu    = yfit.mean(0)
sig   = 1.0*yfit.std(0)
# fit   = abp[0]*x+abp[1]
fit   = abp[0]*xdata

m  = round(abp[0],10)
# b  = round(abp[1],10)
ea = round(abper[0],10)
# eb = round(abper[1],10)

# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

plt.title('$N(OH)_{65}$ and $N(OH)_{67}$', fontsize=30)
plt.ylabel('$N(OH)_{67}$', fontsize=35)
plt.xlabel('$N(OH)_{65}$', fontsize=35)
# plt.xlim(-1.0, 4.0)

plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.text(0.0,4.0, '$Fit: f = ['+str(m)+'\pm'+str(ea) +']\cdot (N_{HI}/10^{20})$', color='blue', fontsize=20)
plt.legend(loc='upper left', fontsize=18)
# for i in range(len(src)):
# 	if(oh[i] > 0):
# 		plt.annotate('('+str(src[i])+', '+ str(xl[i])+', '+str(xb[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
# 	            xytext=(-50.,30.), textcoords='offset points',
# 	            arrowprops=dict(arrowstyle="->"),fontsize=12,
# 	            )
plt.show()
########### END - ODR ############

## ********* Compare N(OH)65 and N(OH)65 from Heiles method and Simple Method ********* ##
xdata = np.array(noh1)
ydata = np.array(xnoh1)

# Error bar for x-axis and y-axis
xerr  = np.array(noh1e)
yerr  = np.array(xnoh1e)

## Fit ##
lguess  = [ 1.]
########### ODR fit ############
# Create a model for Orthogonal distance regression (ODR) fitting.
lin_model = Model(lin_fc)
# Create a RealData object using our initiated data from above.
data      = RealData(xdata, ydata, sx=xerr, sy=yerr)
# Set up ODR with the model and data.
odr       = ODR(data, lin_model, beta0=lguess)
# Run the regression.
out       = odr.run()

print '********* ODR Results *********'
abp   = out.beta
abper = out.sd_beta
for i in range(1):
	print "%s = %03.8f +/- %03.8f" % ('slope',abp[i],abper[i])
## Plot ##
a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
# b     = np.array([ abp[1]-abper[1], abp[1]+abper[1] ])
xfit  = np.linspace(xdata.min(), xdata.max(), 20)
# yfit  = a[:, None] * xfit + b[:, None]
yfit  = a[:, None] * xfit
mu    = yfit.mean(0)
sig   = 1.0*yfit.std(0)
# fit   = abp[0]*x+abp[1]
fit   = abp[0]*xdata

m  = round(abp[0],10)
# b  = round(abp[1],10)
ea = round(abper[0],10)
# eb = round(abper[1],10)

# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

plt.title('$N(OH)_{65-HT}$ and $N(OH)_{65-CB}$', fontsize=30)
plt.ylabel('$N(OH)_{65-HT}$', fontsize=35)
plt.xlabel('$N(OH)_{65-CB}$', fontsize=35)
# plt.xlim(-1.0, 4.0)

plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.text(0.0,4.0, '$Fit: f = ['+str(m)+'\pm'+str(ea) +']\cdot (N_{HI}/10^{20})$', color='blue', fontsize=20)
plt.legend(loc='upper left', fontsize=18)
# for i in range(len(src)):
# 	if(oh[i] > 0):
# 		plt.annotate('('+str(src[i])+', '+ str(xl[i])+', '+str(xb[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
# 	            xytext=(-50.,30.), textcoords='offset points',
# 	            arrowprops=dict(arrowstyle="->"),fontsize=12,
# 	            )
plt.show()
########### END - ODR ############

## ********* Compare N(OH)67 and N(OH)67 from Heiles method and Simple Method ********* ##
xdata = np.array(noh2)
ydata = np.array(xnoh2)

# Error bar for x-axis and y-axis
xerr  = np.array(noh2e)
yerr  = np.array(xnoh2e)

## Fit ##
lguess  = [ 1.]
########### ODR fit ############
# Create a model for Orthogonal distance regression (ODR) fitting.
lin_model = Model(lin_fc)
# Create a RealData object using our initiated data from above.
data      = RealData(xdata, ydata, sx=xerr, sy=yerr)
# Set up ODR with the model and data.
odr       = ODR(data, lin_model, beta0=lguess)
# Run the regression.
out       = odr.run()

print '********* ODR Results *********'
abp   = out.beta
abper = out.sd_beta
for i in range(1):
	print "%s = %03.8f +/- %03.8f" % ('slope',abp[i],abper[i])
## Plot ##
a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
# b     = np.array([ abp[1]-abper[1], abp[1]+abper[1] ])
xfit  = np.linspace(xdata.min(), xdata.max(), 20)
# yfit  = a[:, None] * xfit + b[:, None]
yfit  = a[:, None] * xfit
mu    = yfit.mean(0)
sig   = 1.0*yfit.std(0)
# fit   = abp[0]*x+abp[1]
fit   = abp[0]*xdata

m  = round(abp[0],10)
# b  = round(abp[1],10)
ea = round(abper[0],10)
# eb = round(abper[1],10)

# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

plt.title('$N(OH)_{67-HT}$ and $N(OH)_{67-CB}$', fontsize=30)
plt.ylabel('$N(OH)_{67-HT}$', fontsize=35)
plt.xlabel('$N(OH)_{67-CB}$', fontsize=35)
# plt.xlim(-1.0, 4.0)

plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.text(0.0,4.0, '$Fit: f = ['+str(m)+'\pm'+str(ea) +']\cdot (N_{HI}/10^{20})$', color='blue', fontsize=20)
plt.legend(loc='upper left', fontsize=18)
# for i in range(len(src)):
# 	if(oh[i] > 0):
# 		plt.annotate('('+str(src[i])+', '+ str(xl[i])+', '+str(xb[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
# 	            xytext=(-50.,30.), textcoords='offset points',
# 	            arrowprops=dict(arrowstyle="->"),fontsize=12,
# 	            )
plt.show()
########### END - ODR ############

sys.exit()