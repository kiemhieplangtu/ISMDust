import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import pylab             as pl
import module            as md

from restore             import restore



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

## 79 MS Sources
ms78sc  = md.read_info_ms_79sc(fname = '../result/nhi_lb_thin_cnm_wnm_78src.txt')
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
comsc = []  ## 16 src
difsc = []  ## 14 src
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

print ''
print '=> 16 common src '
for i in range(len(comsc)):
	sc = comsc[i]
	print sc, msdat[sc]['l'], '\t', msdat[sc]['b'], '\t', msdat[sc]['nhi'], '\t', msdat[sc]['nhi_er'], '\t', spdat[sc]['nhi'], '\t', spdat[sc]['nhi_er']

print ''
print 'NHI'
for i in range(len(comsc)):
	sc = comsc[i]
	print sc, spdat[sc]['nhi'], msdat[sc]['nhi'], 100.*(-spdat[sc]['nhi'] + msdat[sc]['nhi'])/msdat[sc]['nhi']

print ''
print 'WNM'
for i in range(len(comsc)):
	sc = comsc[i]
	print sc, spdat[sc]['wnm'], msdat[sc]['wnm'], 100.*(-spdat[sc]['wnm'] + msdat[sc]['wnm'])/msdat[sc]['wnm']

print ''
print 'CNM'
for i in range(len(comsc)):
	sc = comsc[i]
	if(msdat[sc]['cnm'] != 0.):
		print sc, spdat[sc]['cnm'], msdat[sc]['cnm'], 100.*(-spdat[sc]['cnm'] + msdat[sc]['cnm'])/msdat[sc]['cnm']
	else:
		print sc, spdat[sc]['cnm'], msdat[sc]['cnm']

spnhi = np.asarray(spnhi)
msnhi = np.asarray(msnhi)



## NHI vs NHI ##
# plt.plot(msnhi, spnhi, 'r.', ms=12, label='data')
plt.errorbar(msnhi, spnhi, xerr=ms_er, yerr=sp_er, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.plot([0,65],[0,65], 'k--', label='$N^{MS}_{HI} = N^{SP}_{HI}$')
plt.title('$N_{HI}$ SPONGE vs MS', fontsize = 35)
plt.ylabel('$N_{HI}\ SPONGE$', fontsize = 35)
plt.xlabel('$N_{HI}\ MS$', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
# plt.text(0.1, 0.,'14 common sources',color='b',fontsize=18)
plt.legend(loc='upper left', fontsize = 18)
# plt.axvline(x=60., lw=4)
for i in range(len(comsc)):
	sl = str(round(xl[i],2) )
	sb = str(round(xb[i],2) )
	plt.annotate('('+str(comsc[i])+') '+sl+', '+sb, xy=(msnhi[i], spnhi[i]), xycoords='data',
           xytext=(-50.,30.), textcoords='offset points',
           arrowprops=dict(arrowstyle="->"),fontsize=12,
           )
plt.grid()
plt.show()



## Histogram N(HI) ##
plt.hist(np.array(spnhi)/np.array(msnhi), alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=40, lw=3)

plt.title('Histogram of total column density ratio $N^{SP}_{HI}$/$N^{MS}_{HI}$ \n along 14 common lines-of-sight', fontsize=22)
plt.ylabel('Number', fontsize=22,fontweight='bold')
plt.xlabel('$N^{SP}_{HI}/N^{MS}_{HI}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=22)
plt.grid(False)

plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(HI) (<3.0e20)', color='blue', fontsize=20)

plt.tight_layout()
plt.show()



## Histogram N(HI) thin ##
plt.hist(np.array(thinsp)/np.array(thinms), alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=40, lw=3)

plt.title('Histogram of optically thin $N^{SP}_{HI}$/$N^{MS}_{HI}$ \n along 14 common lines-of-sight', fontsize=22)
plt.ylabel('Number', fontsize=22,fontweight='bold')
plt.xlabel('$N^{*SP}_{HI}/N^{*MS}_{HI}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=22)
plt.grid(False)

plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(HI) (<3.0e20)', color='blue', fontsize=20)

plt.tight_layout()
plt.show()



## Histogram CNM ##
plt.hist(np.array(spcnm)/np.array(mscnm), alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=40, lw=3)

plt.title('', fontsize=22)
plt.ylabel('Number', fontsize=22,fontweight='bold')
plt.xlabel('$CNM^{SP}_{HI}/CNM^{MS}_{HI}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=22)
plt.grid(False)

plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(HI) (<3.0e20)', color='blue', fontsize=20)

plt.tight_layout()
plt.show()




## Histogram WNM ##
plt.hist(np.array(spwnm)/np.array(mswnm), alpha=0.5, label='', color='b', ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,2.0), bins=40, lw=3)

plt.title('', fontsize=22)
plt.ylabel('Number', fontsize=22,fontweight='bold')
plt.xlabel('$WNM^{SP}_{HI}/WNM^{MS}_{HI}$', fontsize=22, fontweight='bold')

plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=22)
plt.grid(False)

plt.text(0.05, 15.0, 'Far from 1.0 are l-o-s with very low N(HI) (<3.0e20)', color='blue', fontsize=20)

plt.tight_layout()
plt.show()


########### MPFIT ############
xdata = np.array(msnhi)
ydata = np.array(spnhi)

# Error bar for x-axis and y-axis
xerr = np.array(ms_er)
yerr = np.array(sp_er)

print 'Lenght.....'
print len(xdata)

## Fit ##

# do_linMPfit(xdata, ydata, xerr, yerr, lguess=[1.0])
xfit, yfit, mu, sig, m, ea = md.do_linMPfit(xdata, ydata, xerr, yerr, lguess=[1.0])

# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label=r'$\tau_{353}\ vs\ N^*_{HI}$')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
plt.plot([0,60],[0,60], 'k--', label='$N^*_{H} = N^{LAB}_{HI}$')

plt.title('$N^{SP}_{HI}$ and $N^{MS}_{HI}$ along 14 common lines-of-sight', fontsize=30)
plt.ylabel('$N_{HI} [10^{20}$ cm$^{-2}]$, from 21-SPONGE', fontsize=35)
plt.xlabel('$N_{HI} [10^{20}$ cm$^{-2}]$, MS', fontsize=35)
# plt.xlim(0.0, 2.0)
# plt.ylim(-1.0, 6.0)
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

########### MPFIT - N(HI) thin ############
xdata = np.array(thinms)
ydata = np.array(thinsp)

# Error bar for x-axis and y-axis
xerr = np.array(thinms_er)
yerr = np.array(thinsp_er)

xfit, yfit, mu, sig, m, ea = md.do_linMPfit(xdata, ydata, xerr, yerr, lguess=[1.0])

# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label=r'$\tau_{353}\ vs\ N^*_{HI}$')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
plt.plot([0,60],[0,60], 'k--', label='$N^*_{H} = N^{LAB}_{HI}$')

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
xdata = np.array(thinms)
ydata = np.array(thinsp)

# Error bar for x-axis and y-axis
xerr = np.array(thinms_er)
yerr = np.array(thinsp_er)

xfit, yfit, mu, sig, m, ea = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=[1.0])

# plt.plot(xdata,ydata, 'ok', ls='None', marker='.', lw=1, label='Factor $f = N_{HI}$/$N^*_{HI}$')
plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)
plt.plot([0,60],[0,60], 'k--', label='$x=y$')

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