import os, sys
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md

from restore             import restore

##================= MAIN ========================##
# Cal tau353 from map or not
readmap = False

## Filename of the maps
pth     = os.getenv("HOME")+'/hdata/dust/'
mapFile = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

# Some constants
fk     = 4.77e-7 ## tau353 = 4.77e-27 N(HI) 
plws   = 8.40e-7 ## tau353 = 8.4e-27 N(H) -> Whole sky
pllow  = 6.60e-7 ## tau353 = 6.6e-27 N(H) -> Low NHI

## Read Infor for 93 src
dat   = md.read_info_93src_csv(fname = '../../../doc/observed_src.csv',asarray=True)
xsc   = dat['src']
xl    = dat['l']
xb    = dat['b']
nhi   = dat['nhi']   ## Already in ascending order
nhier = dat['nhier']
xOK   = dat['ok']
tsg67 = dat['tsig67']

fltr  = np.extract([xOK == 1], xOK)
xsc   = np.extract([xOK == 1], xsc)
xl    = np.extract([xOK == 1], xl)
xb    = np.extract([xOK == 1], xb)
nhi   = np.extract([xOK == 1], nhi)
nhier = np.extract([xOK == 1], nhier)
tsg67 = np.extract([xOK == 1], tsg67)
n     = len(fltr)

# to dict. of Infor
infor        = {}
infor['src'] = xsc
infor['l']   = xl
infor['b']   = xb


# tau353 map, err_tau353 map
if(readmap):
	tauMap  = hp.read_map(mapFile, field = 0)
	errMap  = hp.read_map(mapFile, field = 1)
	tau353,\
	tauer   = md.get_tau(tauMap, errMap, infor)
else:
	# Just to run faster
	tau353  = 1e-7*np.array([  21.64003945, 215.25840566,  216.25850422,   92.34412573,  124.53959243, 34.76551456,   39.77825145,    7.84869599,    6.95178528,   65.48004876,
	  207.51052944,   37.28497632,   33.79041118,  137.19225535,   13.04412081,
	   15.96158882,   56.55446785,   13.12147219,  292.5939225,   103.90457646,
	   14.99928999,   16.36758725,   31.24495606,   33.81683655,    8.06284277,
	   17.25587708,   23.24275101,   22.24267519,   57.04321666,    8.13896804,
	  117.95742466,   24.05234454,   23.43321512,   12.96034043])

	tauer   = 1e-8*np.array([   7.09781034,   72.20351108,   81.13784133,   37.0095222,    62.80488378,
	   15.01593943,    9.30117494,    3.15375814,    2.9936075,    34.52255442,
	   62.43462849,   11.19585704,    7.85041721,   40.64422683,    8.53127631,
	    3.11548547,    7.93638719,    9.03687578,  107.52082744,   43.68748705,
	    4.43835688,    7.92562105,    4.17630872,    6.86584372,    2.22341896,
	   27.8553216,     6.29738679,    5.91915068,   14.62017281,    4.90225602,
	   38.90791049,    3.12531938,    1.51442237,    7.29517069])

# Upper limit for N(OH1667)
nohtest = 3.0*tsg67*3.5*2.36634 ## N(OH) = 2.36634*(3*sigma)*3.5(K)*1(Width:km/s)
nhtest  = nhi+20.0*nohtest
sigtest = tau353/nhtest

#############################
#         To fit            #
#############################
xdata = nhi
ydata = tau353
xerr  = nhier
yerr  = tauer

### Correlation: Radiance and N(H) ##
coxy  = md.cov_xy(ydata,xdata)
varx  = md.var(xdata)
vary  = md.var(ydata)
rho   = coxy/varx/vary

print ''
print '********* Pearson Coefficient *********'
print 'Pearson Coeff', md.pearson_coeff(xdata, ydata), ', ', rho
print 'Number of datapoints: ', n
print ''

plt.plot(np.array(xdata)/np.array(varx), np.array(ydata)/np.array(vary), 'r*')
plt.xlabel('X/varx')
plt.ylabel('Y/vary')
plt.show()
### End - Correlation: Tau and NH ##


## Do the Fit and Plot ##
# 1 parameter
lguess = [ 6.8e-7 ]
xfit, yfit, mu, sig, m, ea = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=lguess, xplot=[1.05, 32.0], plot=False)

# 2 params
# lguess = [ 6.8e-7, -4.0e-7]
# xfit, yfit, mu, sig, m, ea, b, eb = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=lguess)


## For plotting
mks    = 8
fts    = 36

plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=mks, markeredgecolor='b', markeredgewidth=1, label='data')
plt.plot([1.0,32.0], [pllow,32.*pllow], 'k-', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='Planck collab. 2014')
plt.plot(xfit, mu, 'b-', mew=2, linewidth=4, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='Linear fit')
plt.plot([1.0,32.0], [fk,32.0*fk], 'r-', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='Fukui et al. 2015')
plt.fill_between(xfit, mu - 2.*sig, mu + 2.*sig, color='0.5', alpha=0.5)

plt.title(r'$\tau_{353}$ and total $N_{H}$', fontsize=30)
plt.ylabel(r'$\tau_{353}$', fontsize=35)
plt.xlabel('$N_{H} [10^{20} cm^{-2}$]', fontsize=35)
plt.xlim(0.8, 90.0)
plt.ylim(5.8e-7, 3.8e-5)

plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
plt.yscale('log')
plt.xscale('log')

plt.text(2.5e-6,1., r'Linear fit: $\tau_{353} = [6.8\pm0.3] \times 10^{-27} \cdot N_{H}$', color='blue', fontsize=20)
plt.text(2.5e-6,1.3, r'Planck collab. 2013: $\tau_{353} = [6.6\pm1.7] \times 10^{-27} \cdot N_{H}$', color='blue', fontsize=20)
plt.text(2.5e-6,1.7, r'Fukui et al. 2015: $\tau_{353} = 4.77 \times 10^{-27} \cdot N_{H}$', color='blue', fontsize=20)
plt.legend(loc='upper left', fontsize=18)
plt.show()
########### END - ODR ############



# Plot for paper Tau353 vs N(H)
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

fig          = plt.figure(figsize=(10,6))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
major_xticks = np.arange(1.e-7, 1.e-4, 1.e-5)
minor_xticks = np.arange(1.e-7, 1.e-4, 5.e-6)
major_yticks = np.arange(1., 80., 10.)                                              
minor_yticks = np.arange(1., 80., 5.0)

ax.errorbar(xdata, ydata*1e6, xerr=xerr, yerr=yerr*1e6, color='k', marker='o', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, capsize=0, label='')
ax.plot(xfit, mu*1e6, 'k-', mew=2, linewidth=2, marker='o', markerfacecolor='b', markersize=0, label=r'$\mathrm{Linear\ fit}$')
ax.plot([1.0,32.0], [pllow*1e6, 32.*pllow*1e6], 'k:', mew=2, linewidth=2, marker='o', markerfacecolor='b', markersize=0, label=r'$\mathrm{PLC2014a}$')
ax.plot([1.0,32.0], [fk*1e6, 32.0*fk*1e6], 'k--', mew=2, linewidth=2, marker='o', markerfacecolor='b', markersize=0, label=r'$\mathrm{Fukui\ et\ al.\ (2015)}$')
ax.fill_between(xfit, mu*1e6 - 3.0*sig*1e6, mu*1e6 + 3.0*sig*1e6, color='0.5', alpha=0.5)

ax.set_ylabel(r'$\mathrm{\tau_{353}\ [10^{-5}]}$', fontsize=26, labelpad=0)
ax.set_xlabel(r'$\mathrm{N_{H}\ [10^{20} cm^{-2}]}$', fontsize=24, labelpad=5)

ax.tick_params(axis='y', labelsize=18)
ax.tick_params(axis='x', labelsize=18)

ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
ax.tick_params(axis='x', labelsize=18, pad=5)
ax.tick_params(axis='y', labelsize=18)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=9)
ax.tick_params(which='minor', length=4)

ax.set_xlim(1.0, 33.0)
ax.set_ylim(5.5e-1, 3.5e1)
ax.set_yscale('log')
ax.set_xscale('log')

# ax.set_xticks([5., 10., 20., 30.])
# ax.set_yticks([1., 5., 10., 20., 30.])
# ax.get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
# ax.get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))

leg = ax.legend(loc='upper left', fontsize=18, numpoints=1)
leg.get_frame().set_linewidth(0.0)

plt.tight_layout()
plt.savefig('tau_vs_nh.eps', bbox_inches='tight', pad_inches=0.1, format='eps', dpi=600)
plt.show()






## PLOT sigma353 vs NHI ##
# From Planck2014a data to overplot
plt.cla()
plt.clf()
plt.close('all')
[nh, sig, sd_sig] = md.read_planck_sigma_vs_nh(fname = 'marc_data/sigma_e353_vs_N_H_Xco1.txt')
nh1    = np.array(nh)
sig1   = np.array(sig)
er1    = np.array(sd_sig)

[nh, sig, sd_sig] = md.read_planck_sigma_vs_nh(fname = 'marc_data/sigma_e353_vs_N_H_Xco2.txt')
nh2    = np.array(nh)
sig2   = np.array(sig)
er2    = np.array(sd_sig)

[nh, sig, sd_sig] = md.read_planck_sigma_vs_nh(fname = 'marc_data/sigma_e353_vs_N_H_Xco3.txt')
nh3    = np.array(nh)
sig3   = np.array(sig)
er3    = np.array(sd_sig)

# This work 
r1     = tau353/nhi
r1er   = md.uncertainty_of_ratio(tau353, nhi, tauer, nhier)

# Sigma353 at low N(HI) from PLC2014a
fltr   = np.extract([(nh2>1.) & (nh2<4.)], sig2)
pl_av  = fltr.mean()

# Sigma353 at low N(HI) from this work
rs     = r1*1e7
rser   = r1er*1e7
fltr   = np.extract([(nhi>1.) & (nhi<4.)], rs)
sg_av  = fltr.mean()

print 'Mean sigma: ', sg_av

# Bin-up data
npoint  = 4
nhiB    = md.bin_up_data(nhi, npoint)
nhierB  = md.bin_up_data(nhier, npoint, error = True)
rsB     = md.bin_up_data(rs, npoint)
rserB   = md.bin_up_data(rser, npoint, error = True)

xx = []
yy = []
for i in range(len(nhiB)):
	if(i==0):
		dx1 = nhiB[i] - nhiB[i+1]
		dx2 = nhiB[i] - nhiB[i+1]
	elif(i==len(nhiB)-1):
		dx1 = nhiB[i-1] - nhiB[i]
		dx2 = nhiB[i-1] - nhiB[i]
	else:
		dx1 = nhiB[i-1] - nhiB[i]
		dx2 = nhiB[i]   - nhiB[i+1]

	xx.append(nhiB[i]+dx1/2.)
	yy.append(rsB[i])

	xx.append(nhiB[i])
	yy.append(rsB[i])

	xx.append(nhiB[i]-dx2/2.)
	yy.append(rsB[i])


## PLOT - sigma353 vs NH ##
mpl.rcParams['axes.linewidth'] = 1.5
fig          = plt.figure(figsize=(16,10))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)

# major_xticks = np.arange(0., 500., 10.)
# minor_xticks = np.arange(0., 500., 10.)
major_yticks = np.arange(0., 22., 2.)                                              
minor_yticks = np.arange(0., 22., 1.)

xerb1, = plt.plot(nh1, sig1, color='b', marker='o', ls='None', markersize=mks, markeredgecolor='b', markeredgewidth=1, label='')
plt.errorbar(nh1, sig1, xerr=nh1*0., yerr=er1*0., color='b', marker='o', ls='None', markersize=mks, markeredgecolor='b', markeredgewidth=1, label='')

xerb2, = plt.plot(nh3, sig3, color='purple', marker='o', ls='None', markersize=mks, markeredgecolor='purple', markeredgewidth=1, label='')
plt.errorbar(nh3, sig3, xerr=nh3*0., yerr=er3*0., color='purple', marker='o', ls='None', markersize=mks, markeredgecolor='purple', markeredgewidth=1, label='')

xerb3, = plt.plot(nh2, sig2, color='0', marker='o', ls='None', markersize=mks, markeredgecolor='k', markeredgewidth=1, label='')
plt.fill_between(nh2, sig2-er2, sig2+er2, color='lightgray', alpha=0.5)

xerb4, = plt.plot(nhi, rs, color='r', marker='d', ls='None', markersize=mks, markeredgecolor='r', markeredgewidth=1, label='')
plt.errorbar(nhi, rs, xerr=nhier, yerr=rser, color='r', marker='d', ls='None', markersize=mks, markeredgecolor='r', markeredgewidth=1, label='')
plt.plot([0.,500.], [sg_av,sg_av], 'k', mew=2, linewidth=2, linestyle='-.', marker='o', markerfacecolor='b', markersize=0, label='')
plt.plot([0.,500.], [pl_av,pl_av], 'k', mew=2, linewidth=2, linestyle='--', marker='o', markerfacecolor='b', markersize=0, label='')

for i in range(n):
	plt.annotate(str(xsc[i]), xy=(nhi[i], rs[i]), xycoords='data',
			            xytext=(-20.,10.), textcoords='offset points',
			            arrowprops=dict(arrowstyle="->"),fontsize=12,
			            )

# Bin-up #
# plt.plot(xx, yy, 'k', mew=2, linewidth=2, linestyle='-', marker='o', markerfacecolor='b', markersize=0, label='')
# plt.errorbar(nhiB, rsB, xerr=nhierB, yerr=rserB, color='r', marker='o', ls='None', markersize=mks, markeredgecolor='r', markeredgewidth=1, label='')

# Test
plt.plot(nhtest, sigtest*1e7, color='violet', marker='d', ls='None', markersize=mks, markeredgecolor='violet', markeredgewidth=1, label='')

plt.title('', fontsize=0)
plt.ylabel(r'$\mathrm{\sigma_{353} = \tau_{353}/N_{H}\ [10^{-27}cm^{2}H^{-1}]} $', fontsize=fts, fontweight='normal')
plt.xlabel(r'$\mathrm{N_{H} [10^{20} cm^{-2}]}$', fontsize=fts, fontweight='normal')

plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=15)
plt.xscale('log')
                                     
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=22, pad=7)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=12)
plt.tick_params(which='minor', length=6)
plt.grid(False)

plt.xlim(0.7, 500.0)
plt.ylim(2.,20.)

axbox = ax.get_position()
leg   = plt.legend([xerb4, xerb1, xerb3, xerb2], [r'$\mathrm{This\ work}$',\
 r'$\mathrm{X_{CO}}$$\mathrm{{=}1}$$\times$$\mathrm{10}$$^{20}$ $\mathrm{(PLC2014a)}$',\
  r'$\mathrm{X_{CO}}$$\mathrm{{=}2}$$\times$$\mathrm{10}$$^{20}$ $\mathrm{(PLC2014a)}$',\
   r'$\mathrm{X_{CO}}$$\mathrm{{=}3}$$\times$$\mathrm{10}$$^{20}$ $\mathrm{(PLC2014a)}$'],\
    fontsize=23, loc=(axbox.x0-0.11, axbox.y0+0.6), numpoints=1, handletextpad=-0.5)
leg.get_frame().set_linewidth(0.0)

plt.tight_layout()

plt.savefig('sig353_vs_nh.eps', bbox_inches='tight', pad_inches=0.1, format='eps', dpi=150)
# plt.savefig('sig353_vs_nh.png', bbox_inches='tight', pad_inches=0.03, format='png', dpi=100)
plt.show()
## END - PLOT ##

sys.exit()