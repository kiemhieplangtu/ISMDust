import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import operator
import copy

from restore             import restore
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

## Linear ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
# def tb_exp(x,y,bg,tau,v0,wid,tex,cont,err=1.):
def myfunc(p, fjac=None, x=None, y=None, err=None):
	# model  = p[0] * x + p[1]
	model  = p[0] * x
	status = 0
	return [status, (y - model) / err]

## Read info of 30 SPONGE sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_info_sponge_30src(fname = '../../oh/sponge/sub_data/30src_claire.txt'):
	cols = ['src','l', 'b', 'cnm','cnm_er','wnm','wnm_er','nhi','nhi_er','thin','thin_er']
	fmt  = ['s',  'f', 'f', 'f',    'f',   'f',  'f',     'f',   'f',     'f',    'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

## Read info of 30 SPONGE sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_info_ms_78src(fname = '../rearrange/nhi_lb_thin_78src.txt'):
	cols = ['idx','src','l', 'b', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	return dat

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
		ret[sc[i]]           = {}
		ret[sc[i]]['l']      = xl[i]
		ret[sc[i]]['b']      = xb[i]
		ret[sc[i]]['nhi']    = hi[i]
		ret[sc[i]]['nhi_er'] = hier[i]
		ret[sc[i]]['thin']   = thin[i]
		ret[sc[i]]['thiner'] = thiner[i]
		ret[sc[i]]['cnm']    = cnm[i]
		ret[sc[i]]['cnm_er'] = cnmer[i]
		ret[sc[i]]['wnm']    = wnm[i]
		ret[sc[i]]['wnm_er'] = wnmer[i]

	return ret

## ======== MAIN ============= ##
## 21SPONGE Claire 30 Sources
sp30sc  = read_info_sponge_30src(fname = '../../oh/sponge/sub_data/30src_claire.txt')
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
ms78sc = read_info_ms_78src(fname = '../result/nhi_lb_thin_cnm_wnm_78src.txt')
sc2     = ms78sc['src']
xl2     = ms78sc['l']
xb2     = ms78sc['b']
hi2     = ms78sc['nhi']
hi2er   = ms78sc['nhi_er']
thin2   = ms78sc['thin']
thin2er = ms78sc['thin_er']
cnm2    = ms78sc['cnm']
cnm2er  = ms78sc['cnm_er']
wnm2    = ms78sc['wnm']
wnm2er  = ms78sc['wnm_er']
msdat   = read2dict(ms78sc)

## Common & different Sources
comsc = []  ## 14 src
difsc = []  ## 16 src
for sc in spdat:
	if (sc in msdat):
		comsc.append(sc)
	else:
		difsc.append(sc)

spnhi     = []
msnhi     = []
sp_er     = []
ms_er     = []
spwnm     = []
mswnm     = []
spwnmer   = []
mswnmer   = []
spcnm     = []
spcnmer   = []
mscnm     = []
mscnmer   = []
xl        = []
xb        = []
thinms    = []
thinms_er = []
thinsp    = []
thinsp_er = []
for i in range(len(comsc)):
	sc = comsc[i]
	# print sc, spdat[sc]['l'], msdat[sc]['l'], spdat[sc]['l'] - msdat[sc]['l']
	print sc, spdat[sc]['nhi'], msdat[sc]['nhi'], 100.*(spdat[sc]['nhi'] - msdat[sc]['nhi'])/msdat[sc]['nhi']
	spnhi.append(spdat[sc]['nhi'])
	msnhi.append(msdat[sc]['nhi'])
	spwnm.append(spdat[sc]['wnm'])
	mswnm.append(msdat[sc]['wnm'])
	spcnm.append(spdat[sc]['cnm'])
	spcnmer.append(spdat[sc]['cnm_er'])
	mscnm.append(msdat[sc]['cnm'])
	mscnmer.append(msdat[sc]['cnm_er'])

	spwnm.append(spdat[sc]['wnm'])
	spwnmer.append(spdat[sc]['wnm_er'])
	mswnm.append(msdat[sc]['wnm'])
	mswnmer.append(msdat[sc]['wnm_er'])

	ms_er.append(msdat[sc]['nhi_er'])
	sp_er.append(spdat[sc]['nhi_er'])
	xl.append(msdat[sc]['l'])
	xb.append(msdat[sc]['b'])
	thinms.append(msdat[sc]['thin'])
	thinms_er.append(msdat[sc]['thiner'])
	thinsp.append(spdat[sc]['thin'])
	thinsp_er.append(spdat[sc]['thiner'])

print 'NHI'
for i in range(len(comsc)):
	sc = comsc[i]
	print sc, spdat[sc]['nhi'], msdat[sc]['nhi'], 100.*(-spdat[sc]['nhi'] + msdat[sc]['nhi'])/msdat[sc]['nhi']

print 'WNM'
for i in range(len(comsc)):
	sc = comsc[i]
	print sc, spdat[sc]['wnm'], msdat[sc]['wnm'], 100.*(-spdat[sc]['wnm'] + msdat[sc]['wnm'])/msdat[sc]['wnm']

print 'CNM'
for i in range(len(comsc)):
	sc = comsc[i]
	if(msdat[sc]['cnm'] != 0.):
		print sc, spdat[sc]['cnm'], msdat[sc]['cnm'], 100.*(-spdat[sc]['cnm'] + msdat[sc]['cnm'])/msdat[sc]['cnm']
	else:
		print sc, spdat[sc]['cnm'], msdat[sc]['cnm']

spnhi = np.asarray(spnhi)
msnhi = np.asarray(msnhi)

########### Low NHI < 3.0e20 ############
lownhi     = []
lownhi_er  = []
low_sp     = []
low_sp_er  = []
for i in range(len(thinms)):
	if(thinms[i] < 5.0):
		lownhi.append(thinms[i])
		lownhi_er.append(thinms_er[i])
		low_sp.append(thinsp[i])
		low_sp_er.append(thinsp_er[i])

########### ODR fit ############
ydata = np.array(low_sp)
xdata = np.array(lownhi)

# Error bar for x-axis and y-axis
yerr = np.array(low_sp_er)
xerr = np.array(lownhi_er)


########### MPFIT - N(HI) thin ############
## Fit ##
lguess  = [1.0]

npar    = len(lguess)
guessp  = np.array(lguess, dtype='float64')
plimd   = [[False,False]]*npar
plims   = [[0.,0.]]*npar
parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]}
pname   = ['slope','offset']
pfix    = [False]*npar

parinfo = []
for i in range(len(guessp)):
	parinfo.append(copy.deepcopy(parbase))

for i in range(len(guessp)):
	parinfo[i]['value']   = guessp[i]
	parinfo[i]['fixed']   = pfix[i]
	parinfo[i]['parname'] = pname[i]
	parinfo[i]['limited'] = plimd[i]

x  = xdata.astype(np.float64)
y  = ydata.astype(np.float64)
er = yerr.astype(np.float64)

fa = {'x':x, 'y':y, 'err':er}
mp = mpfit(myfunc, guessp, parinfo=parinfo, functkw=fa, quiet=True)

## ********* Results ********* ##
print '********* Results *********'
abp   = mp.params
abper = mp.perror
for i in range(len(parinfo)):
	print "%s = %03.8f +/- %03.8f" % (parinfo[i]['parname'],abp[i],abper[i])
## Plot ##
a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
# b     = np.array([ abp[1]-abper[1], abp[1]+abper[1] ])
xfit  = np.linspace(xdata.min(), xdata.max(), 20)
# yfit  = a[:, None] * xfit + b[:, None]
yfit  = a[:, None] * xfit
mu    = yfit.mean(0)
sig   = 1.0*yfit.std(0)
# fit   = abp[0]*x+abp[1]
fit   = abp[0]*x

m  = round(abp[0],10)
ea = round(abper[0],10)
# b  = round(abp[1],10)
# eb = round(abper[1],10)

# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label=r'$\tau_{353}\ vs\ N^*_{HI}$')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
# plt.plot([0,60],[0,60], 'k--', label='$N^*_{H} = N^{LAB}_{HI}$')

plt.title('$N^{*SP}_{HI}$ and $N^{*MS}_{HI}$ along 14 common lines-of-sight', fontsize=30)
plt.ylabel('$N^{thin}_{HI} [10^{20}$ cm$^{-2}]$, from 21-SPONGE', fontsize=35)
plt.xlabel('$N^{thin}_{HI} [10^{20}$ cm$^{-2}]$, MS', fontsize=35)
plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.text(0.0, 40.0, '$f = ['+str(m)+'\pm'+str(ea) +']\cdot (N^*_{HI}/10^{20})$', color='blue', fontsize=20)
plt.legend(loc='upper left', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')
# for i in range(len(src)):
# 	plt.annotate('('+str(src[i])+', '+ str(xl[i])+', '+str(xb[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
#             xytext=(-50.,30.), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->"),fontsize=12,
#             )
plt.show()
########### END - MPFIT ############

########### ODR fit ############
# Create a model for Orthogonal distance regression (ODR) fitting.
lin_model = Model(lin_fc)
# Create a RealData object using our initiated data from above.
data      = RealData(xdata, ydata, sx=xerr, sy=yerr)
# Set up ODR with the model and data.
odr       = ODR(data, lin_model, beta0=[1.])
# Run the regression.
out       = odr.run()

## ********* Results ********* ##
print '********* ODR Results *********'
abp   = out.beta
abper = out.sd_beta
for i in range(len(parinfo)):
	print "%s = %03.8f +/- %03.8f" % (parinfo[i]['parname'],abp[i],abper[i])
## Plot ##
a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
# b     = np.array([ abp[1]-abper[1], abp[1]+abper[1] ])
xfit  = np.linspace(xdata.min(), xdata.max(), 20)
# yfit  = a[:, None] * xfit + b[:, None]
yfit  = a[:, None] * xfit
mu    = yfit.mean(0)
sig   = 1.0*yfit.std(0)
# fit   = abp[0]*x+abp[1]
fit   = abp[0]*x

m  = round(abp[0],10)
# b  = round(abp[1],10)
ea = round(abper[0],10)
# eb = round(abper[1],10)

# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
# plt.plot([0,60],[0,60], 'k--', label='$x=y$')

plt.title('$N^{*SP}_{HI}$ and $N^{*MS}_{HI}$ along 14 common lines-of-sight', fontsize=30)
plt.ylabel('$N^{thin}_{HI} [10^{20}$ cm$^{-2}]$, from 21-SPONGE', fontsize=35)
plt.xlabel('$N^{thin}_{HI} [10^{20}$ cm$^{-2}]$, MS', fontsize=35)
# plt.xlim(0.0, 2.0)
# plt.ylim(-1.0, 6.0)
plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
# plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0))

plt.text(0.0, 40.0, '$y = ['+str(m)+'\pm'+str(ea) +']\cdot x$', color='blue', fontsize=20)
plt.legend(loc='upper left', fontsize=18)
# plt.savefig("test.png",bbox_inches='tight')
# for i in range(len(src)):
# 	plt.annotate('('+str(src[i])+', '+ str(xl[i])+', '+str(xb[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
#             xytext=(-50.,30.), textcoords='offset points',
#             arrowprops=dict(arrowstyle="->"),fontsize=12,
#             )
plt.show()
########### END - ODR ############

sys.exit()

## Plot DIFFERENCE ##
y = 100.*(msnhi-spnhi)/msnhi
plt.plot(y, 'r*-', ms=12, label='$[(N^{MS}_{HI}-N^{SP}_{HI})/N^{MS}_{HI}] (\%)$')
plt.grid()
plt.title('$N_{HI}$ Difference, SPONGE vs MS', fontsize = 35)
plt.ylabel('$[(N^{MS}_{HI}-N^{SP}_{HI})/N^{MS}_{HI}] (\%)$', fontsize = 35)
plt.xlabel('Index', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(0.1, 0.,'14 common sources',color='b',fontsize=18)
plt.legend(loc='upper left', fontsize = 18)
# plt.axvline(x=60., lw=4)
for i in range(len(comsc)):
	sl = str(round(xl[i],2) )
	sb = str(round(xb[i],2) )
	plt.annotate('('+str(comsc[i])+') '+sl+', '+sb, xy=(i, y[i]), xycoords='data',
           xytext=(-50.,30.), textcoords='offset points',
           arrowprops=dict(arrowstyle="->"),fontsize=18,
           )
plt.show()

## WNM
plt.subplot(2,1,1)
spwnm = np.asarray(spwnm)
mswnm = np.asarray(mswnm)
y     = 100.*(mswnm-spwnm)/mswnm
plt.plot(y, 'r*-', ms=12, label='$[(N^{MS}_{wnm}-N^{SP}_{wnm})/N^{MS}_{wnm}] (\%)$')
plt.grid()
plt.title('$N_{wnm}$ Difference, SPONGE vs MS', fontsize = 35)
plt.ylabel('$[(N^{MS}_{wnm}-N^{SP}_{wnm})/N^{MS}_{wnm}] (\%)$', fontsize = 35)
# plt.xlabel('Index', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(0.1, 0.,'14 common sources',color='b',fontsize=18)
plt.legend(loc='upper left', fontsize = 18)
# plt.axvline(x=60., lw=4)
for i in range(len(comsc)):
	sl = str(round(xl[i],2) )
	sb = str(round(xb[i],2) )
	plt.annotate('('+str(comsc[i])+') '+sl+', '+sb, xy=(i, y[i]), xycoords='data',
           xytext=(-50.,-20.), textcoords='offset points',
           arrowprops=dict(arrowstyle="->"),fontsize=12,
           )

## CNM
plt.subplot(2,1,2)
spcnm = np.asarray(spcnm)
mscnm = np.asarray(mscnm)
y     = 100.*(mscnm-spcnm)/mscnm
plt.plot(y, 'r*-', ms=12, label='$[(N^{MS}_{cnm}-N^{SP}_{cnm})/N^{MS}_{cnm}] (\%)$')
plt.grid()
plt.title('$N_{cnm}$ Difference, SPONGE vs MS', fontsize = 35)
plt.ylabel('$[(N^{MS}_{cnm}-N^{SP}_{cnm})/N^{MS}_{cnm}] (\%)$', fontsize = 35)
plt.xlabel('Index', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.text(0.1, 0.,'14 common sources',color='b',fontsize=18)
plt.legend(loc='upper left', fontsize = 18)
# plt.axvline(x=60., lw=4)
for i in range(len(comsc)):
	sl = str(round(xl[i],2) )
	sb = str(round(xb[i],2) )
	plt.annotate('('+str(comsc[i])+') '+sl+', '+sb, xy=(i, y[i]), xycoords='data',
           xytext=(-50.,-20.), textcoords='offset points',
           arrowprops=dict(arrowstyle="->"),fontsize=12,
           )
plt.show()