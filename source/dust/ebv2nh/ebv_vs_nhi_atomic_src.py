import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np
import healpy            as hp
import module            as md

from   restore           import restore

#================= MAIN ========================#
# Cal EBV from map or not
readmap = False

# data file
pth     = os.getenv("HOME")+'/hdata/dust/'
mapFile = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'  ## E(B-V) from SFD et al. 1998, IRAS ~6'

## Infor for 93 src
dat   = md.read_info_93src_csv(fname = '../../../doc/observed_src.csv',asarray=True)
xsc   = dat['src']
xl    = dat['l']
xb    = dat['b']
nhi   = dat['nhi']
nhier = dat['nhier']
ebv   = dat['ebv']
ebver = dat['ebver']
xOK   = dat['ok']

fltr  = np.extract([xOK == 1], xOK)
xsc   = np.extract([xOK == 1], xsc)
xl    = np.extract([xOK == 1], xl)
xb    = np.extract([xOK == 1], xb)
nhi   = np.extract([xOK == 1], nhi)
nhier = np.extract([xOK == 1], nhier)
ebv   = np.extract([xOK == 1], ebv)
ebver = np.extract([xOK == 1], ebver)
n     = len(fltr)

# to dict. of Infor
infor        = {}
infor['src'] = xsc
infor['l']   = xl
infor['b']   = xb

########### FIT ############
ydata = ebv
xdata = nhi

# Error bar for x-axis and y-axis
yerr = ebver
xerr = nhier

print xdata.shape
print ydata.shape

print ''
### Correlation: Radiance and NH ##
coxy      = md.cov_xy(ydata,xdata)
varx      = md.var(xdata)
vary      = md.var(ydata)
rho       = coxy/varx/vary

print ''
print '********* Pearson Coefficient *********'
print 'Pearson Coeff', md.pearson_coeff(xdata, ydata), ', ', rho
print 'Number of datapoints: ', n
print ''

plt.plot(xdata/varx, ydata/vary, 'r*')
plt.xlabel('X/varx')
plt.ylabel('Y/vary')
plt.show()
### End - Correlation: Radiance and NH ##

## Fit ODR for data from Schlafly 2011 ##	
# xfit, yfit, mu, sig, m, ea, b, eb = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=[0.017, 0.01], xplot=[1.05, 32.], plot=False)
xfit, yfit, mu, sig, m, ea = md.do_linODRfit(ydata, xdata, yerr, xerr, lguess=[0.017], xplot=[1.05, 32.], plot=False)

sys.exit()


### E(B-V) from Planck R1.2 ###
# TFORM1  = 'EBV       '          /data format of field: 4-byte REAL
# TTYPE2  = 'N_OBS   '           /label for field   2
if(readmap):
	tau_map    = hp.read_map(mapFile, field = 0)
	tauer_map  = hp.read_map(mapFile, field = 1)
	ci_map     = hp.read_map(mapFile, field = 2)
	ci, cier   = md.cal_ebv_from_planckR12(tau_map, tauer_map, ci_map, infor)
else:
	# just to run quickly
	ci   = np.array([0.0331172943115,  0.38322353363,  0.263280242682,  0.111670605838,  0.162273988128,  0.0510572008789,  0.0647942721844,  0.0217490531504,  0.015030387789,  0.0986406430602,  0.273433506489,  0.0521479435265,  0.0660882815719,  0.19742205739,  0.0281760767102,  0.0313837938011,  0.081213131547,  0.0246310308576,  0.274701237679,  0.138853147626,  0.0351148806512,  0.0441353805363,  0.0402345992625,  0.0474565215409,  0.0209110565484,  0.0340617895126,  0.0406733937562,  0.031562063843,  0.109481349587,  0.01621058397,  0.216225162148,  0.0405353605747,  0.0365297757089,  0.0335537865758])
	cier = np.array([0.00127455984039,  0.0149923237275,  0.0112104715642,  0.00500854760775,  0.00881154778939,  0.00243309959638,  0.00199933268341,  0.000977492381742,  0.00071449850667,  0.00556688385465,  0.00989905520171,  0.00188531638339,  0.00203176536055,  0.0070716566423,  0.00192815166145,  0.000880070382876,  0.00199314486667,  0.00176736662491,  0.0115104826714,  0.0064730510139,  0.00125678982016,  0.00231453825787,  0.000972351593851,  0.00135695807262,  0.000713992210569,  0.00554103336749,  0.00137297317176,  0.00105323175914,  0.00356829215426,  0.0010295023075,  0.00835585910832,  0.000971349898213,  0.00077245910278,  0.00200588039363])

# data
xx   = nhi
yy   = ci

xxer = nhier
yyer = cier

########### FIT ############
xdat = xx
ydat = yy

# Error bar for x-axis and y-axis
xer  = xxer
yer  = yyer

print ''
### Correlation: Radiance and NH ##
coxy = md.cov_xy(ydat,xdat)
varx = md.var(xdat)
vary = md.var(ydat)
rho  = coxy/varx/vary

print ''
print '********* Pearson Coefficient *********'
print 'Pearson Coeff', md.pearson_coeff(xdat, ydat), ', ', rho
print 'Number of datapoints: ', n
print ''

plt.plot(np.array(xdat)/np.array(varx), np.array(ydat)/np.array(vary), 'r*')
plt.xlabel('X/varx')
plt.ylabel('Y/vary')
plt.show()

## Fit ODR for data from Placnk 2014 R1.2 ##	
# xFit, yFit, Mu, Sig, M, Ea, B, Eb = md.do_linODRfit(xdat, ydat, xer, yer, lguess=[0.017, 0.01], xplot=[1.05, 32.], plot=False)
xFit, yFit, Mu, Sig, M, Ea = md.do_linODRfit(xdat, ydat, xer, yer, lguess=[0.017], xplot=[1.05, 32.], plot=False)
### End - Correlation: Radiance and NH ###
### END - E(B-V) from Planck R1.2      ###


### PLOT FIGURES ###
mks = 8
fts = 36

mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

f,(ax1, ax2) = plt.subplots(2, 1, gridspec_kw = {'wspace':0, 'hspace':0.02}, figsize=(10,8))

major_xticks = np.arange(0., 0.5, 0.1)
minor_xticks = np.arange(0., 0.5, 0.05)
major_yticks = np.arange(0., 80., 10.)                                              
minor_yticks = np.arange(0., 80., 5.0)

ax1.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='k', marker='o', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, capsize=0, label='data')
ax1.plot(xfit, mu, '-k', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='k', markersize=0, label='ODR linear fit')
ax1.fill_between(xfit, mu-3.0*sig, mu+ 3.0*sig, color='0.5', alpha=0.5)
xerb1, = ax1.plot(xdata, ydata, color='k', marker='o', ls='None', markersize=5, markeredgecolor='k', markeredgewidth=1, label='') 

ax1.set_ylabel(r'$\mathrm{E(B{-}V)\ [mag]}$', fontsize=18)

ax1.set_xticks(major_xticks)                                                       
ax1.set_xticks(minor_xticks, minor=True)                                           
ax1.set_yticks(major_yticks)                                                       
ax1.set_yticks(minor_yticks, minor=True)
ax1.tick_params(axis='x', pad=0, labelbottom='off')
ax1.tick_params(axis='y', labelsize=18)
ax1.tick_params(which='both', width=1.5)
ax1.tick_params(which='major', length=8)
ax1.tick_params(which='minor', length=4)

ax1.set_xlim(1.0, 33.)
ax1.set_ylim(0.007, 0.5)

ax1.set_xscale('log')
ax1.set_yscale('log')

# ax1.set_xticks([0.01, 0.05, 0.1, 0.3])
# ax1.set_yticks([1, 5, 10, 30])
# ax1.get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
# ax1.get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))

axbox = ax1.get_position()
leg   = ax1.legend([xerb1], [r'$\mathrm{E(B{-}V)\ from\ Schlafly\ et\ al.\ (2011)}$'], \
	 fontsize=13, loc=(axbox.x0-0.12, axbox.y0+0.33), numpoints=1, handletextpad=0.)
leg.get_frame().set_linewidth(0.0)


### Ax2 ###
ax2.errorbar(xdat, ydat, xerr=xer, yerr=yer, color='k', marker='h', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, capsize=0, label='data')
ax2.plot(xFit, Mu, '-k', mew=2, linewidth=2, linestyle='solid', marker='h', markerfacecolor='k', markersize=0, label='ODR linear fit')
ax2.fill_between(xFit, Mu-3.0*Sig, Mu+3.0*Sig, color='0.5', alpha=0.5)
xerb1, = ax2.plot(xdat, ydat, color='k', marker='h', ls='None', markersize=5, markeredgecolor='k', markeredgewidth=1, label='') 

# plt.title('', fontsize=30)
ax2.set_ylabel(r'$\mathrm{E(B{-}V)\ [mag]}$', fontsize=18)
ax2.set_xlabel(r'$\mathrm{N_{H}\ [10^{20} cm^{-2}]}$', fontsize=20)

ax2.set_xticks(major_xticks)                                                       
ax2.set_xticks(minor_xticks, minor=True)                                           
ax2.set_yticks(major_yticks)                                                       
ax2.set_yticks(minor_yticks, minor=True)
ax2.tick_params(axis='x', labelsize=18, pad=5)
ax2.tick_params(axis='y', labelsize=18)
ax2.tick_params(which='both', width=1.5)
ax2.tick_params(which='major', length=8)
ax2.tick_params(which='minor', length=4)

ax2.set_xlim(1.0, 33.)
ax2.set_ylim(0.007, 0.5)

ax2.set_xscale('log')
ax2.set_yscale('log')

# ax2.set_xticks([0.01, 0.05, 0.1, 0.3])
# ax2.set_yticks([1, 5, 10, 30])
# ax2.get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.2f'))
# ax2.get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))

axbox = ax2.get_position()
leg   = ax2.legend([xerb1], [r'$\mathrm{E(B{-}V)\ from\ PLC2014a}$'], \
	 fontsize=13, loc=(axbox.x0-0.12, axbox.y0+0.73), numpoints=1, handletextpad=0.)
leg.get_frame().set_linewidth(0.0)

plt.savefig('ebv_vs_nh.eps', bbox_inches='tight', pad_inches=0.08, format='eps', dpi=600)
plt.show()
########### END - ODR ############



# for i in range(n):
# 	nh = m*xdata[i]+b
# 	print i, xsc[i],'     \t', xdata[i],'\t', ydata[i],'\t', nh,'\t', nh-ydata[i]  ## compare/nh_from_ebv2011.txt