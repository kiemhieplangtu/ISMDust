import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md

from numpy               import array
from restore             import restore
from plotting            import cplot
from mpfit               import mpfit


##================= MAIN ========================##
# Cal tau353 from map or not
readmap = False

## Filename of the map
pth     = os.getenv("HOME")+'/hdata/dust/'
mapFile = pth + 'HFI_CompMap_ThermalDustModel_2048_R1.20.fits'

## Infor for 93 src
dat   = md.read_info_93src_csv(fname = '../../../doc/observed_src.csv',asarray=True)
xsc   = dat['src']
xl    = dat['l']
xb    = dat['b']
nhi   = dat['nhi']
nhier = dat['nhier']
xOK   = dat['ok']

fltr  = np.extract([xOK == 1], xOK)
xsc   = np.extract([xOK == 1], xsc)
xl    = np.extract([xOK == 1], xl)
xb    = np.extract([xOK == 1], xb)
nhi   = np.extract([xOK == 1], nhi)
nhier = np.extract([xOK == 1], nhier)
n     = len(fltr)

# to dict. of Infor
infor        = {}
infor['src'] = xsc
infor['l']   = xl
infor['b']   = xb


# Radiance map, err_R map and resolution #
if(readmap):
	tauMap    = hp.read_map(mapFile, field = 0)
	tauErrMap = hp.read_map(mapFile, field = 1)
	rMap      = hp.read_map(mapFile, field = 3)
	r, rer    = md.get_radiance(tauMap, tauErrMap, rMap, infor)
else:
	# Just to run 
	r   = 1e-8*np.array([6.12065420569,  70.9746871053,  48.7563227125,  20.6769968258,  30.0596269653,  9.46087368447,  12.0033092799,  4.02382482889,  2.78365739348,  18.2798800097,  50.6284209223,  9.63273407706,  12.2375382716,  36.5614113207,  5.10232602835,  5.80633177094,  15.0493121964,  4.56085800238,  50.8637242547,  25.7145359228,  6.50484039966,  8.17065952674,  7.44642463246,  8.78252137682,  3.85327751928,  6.34381294162,  7.53960023303,  5.84083359456,  20.2823670747,  3.00155136301,  40.045466676,  7.51680673261,  6.76277664979,  6.21246840637])
	rer = 1e-9*np.array([2.56508149747,  30.159172782,  22.2800827312,  9.88776020327,  17.0670905463,  4.77383604987,  4.20505508257,  1.92770920741,  1.40151053259,  10.7528041501,  20.1612009171,  3.83159831076,  4.27489039012,  14.4323318649,  3.5927447391,  1.89174771968,  4.45789812588,  3.35888844065,  22.9220235171,  12.7238286096,  2.56601222703,  4.49407313455,  2.18268981974,  2.90320784179,  1.4627091877,  10.3733662084,  2.83574494678,  2.17661082351,  7.41743072625,  1.97017323985,  16.8406313173,  2.19069303902,  1.81753445554,  3.85419226002])


#############################
#         To fit            #
#############################
xdata     = nhi
ydata     = r
xerr      = nhier
yerr      = rer

### Correlation: Radiance and NH ##
coxy      = md.cov_xy(xdata,ydata)
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

########### MPFIT ############
ydata  = 1.e-4*ydata*1.0e11   ## radiance W/m2/sr -> w/cm2/sr, xdata*1e11 => 1e7
xdata  = xdata                ## 1e20

# Error bar for x-axis and y-axis
yerr   = 1.e-4*yerr*1.0e11     ## radiance W/m2/sr -> w/cm2/sr, xdata*1e11 => 1e7
xerr   = xerr                  ## 1e20

lguess = [0.2]
########### MPFIT fit ############
# xfit, yfit, mu, sig, m, ea, b, eb = md.do_linMPfit(xdata, ydata, xerr, yerr, lguess=[3.0, 1.0])

########### ODR fit ############
# xfit, yfit, mu, sig, m, ea = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=[0.2])
xfit, yfit, mu, sig, m, ea = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=lguess, xplot=[1.05, 32.0], plot=False)

print 'Fit Results:'
print m
print ea


## For plotting
mks    = 8
fts    = 36

plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=mks, markeredgecolor='b', markeredgewidth=1, label='data')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.plot(xdata, xdata*0.21, 'r-', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu-3.0*sig, mu+3.0*sig, color='0.5', alpha=0.5)

plt.ylabel('$Radiance\ [10^{-11} Wcm^{-2}sr^{-1}]$', fontsize=35)
plt.xlabel('$N_{HI} [10^{20} cm^{-2}]$', fontsize=35)

# plt.xlim(0.2, 10.0)
# plt.ylim(-1.0, 6.0)

plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
plt.yscale('log')
plt.xscale('log')
plt.legend(loc='upper left', fontsize=18)
plt.show()
########### END - ODR ############

for i in range(n):
	nh = 4.046e31*xdata[i]/1e11+0.09269e20
	nh = nh/1e20
	print i, xsc[i],'     \t', xdata[i],'\t', ydata[i],'\t', nh,'\t', nh-ydata[i]


# Plot for paper R vs N(H)
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

fig          = plt.figure(figsize=(10,6))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
major_xticks = np.arange(0., 12., 10.)
minor_xticks = np.arange(0., 12., 5.)
major_yticks = np.arange(0., 80., 10.)                                              
minor_yticks = np.arange(0., 80., 5.0)

plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='k', marker='o', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, capsize=0, label='data')
plt.plot(xfit, mu, '-k', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='k', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu-3.0*sig, mu+3.0*sig, color='0.5', alpha=0.5)

plt.xlabel(r'$\mathrm{N_{H}\ [10^{20} cm^{-2}]}$', fontsize=24)
plt.ylabel(r'$\mathrm{\mathcal{R}\ [10^{-11} Wcm^{-2}sr^{-1}]}$', fontsize=24)

ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
ax.tick_params(axis='x', labelsize=18, pad=5)
ax.tick_params(axis='y', labelsize=18)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=9)
ax.tick_params(which='minor', length=4)

ax.set_xlim(1.0, 33.)
ax.set_ylim(0.2, 10.)

ax.set_xscale('log')
ax.set_yscale('log')

# ax.set_xticks([5., 10., 20., 30.])
# ax.set_yticks([0.5, 1., 5., 10.])
# ax.get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
# ax.get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))

plt.tight_layout()
plt.savefig('r_vs_nh.eps', bbox_inches='tight', pad_inches=0.08, format='eps', dpi=600)
plt.show()








## PLOT sigma vs NHI ##
# From Planck data to overplot
[nh, sig, sd_sig] = md.read_planck_sigma_vs_nh(fname = 'marc_data/L_H_vs_N_H_Xco1.txt')
nh1     = np.array(nh)
sig1    = np.array(sig)
er1     = np.array(sd_sig)

[nh, sig, sd_sig] = md.read_planck_sigma_vs_nh(fname = 'marc_data/L_H_vs_N_H_Xco2.txt')
nh2     = np.array(nh)
sig2    = np.array(sig)
er2     = np.array(sd_sig)

[nh, sig, sd_sig] = md.read_planck_sigma_vs_nh(fname = 'marc_data/L_H_vs_N_H_Xco3.txt')
nh3     = np.array(nh)
sig3    = np.array(sig)
er3     = np.array(sd_sig)

# This work
r1      = 4.*3.14159*r/nhi
r1er    = md.uncertainty_of_ratio(r, nhi, rer, nhier)

# CNM fraction ??
# fcnm1   = cnm/nhi
# f1er    = md.uncertainty_of_ratio(cnm, nhi, cnmer, nhier)

## Mean L at low N(HI)
fltr    = np.extract([(nh2>1.) & (nh2<4.)], sig2)
pl_av   = fltr.mean()

print pl_av

## Mean L at low N(HI) from this work
r1      = r1*1e7
r1er    = r1er*1e7*4.*3.14159
nh      = nhi

fltr    = np.extract([(nh>1.) & (nh<4.)], r1)
r_av    = fltr.mean()

print 'r_av', r_av


# Bin-up data
npoint  = 4
nhiB    = md.bin_up_data(nhi, npoint)
nhierB  = md.bin_up_data(nhier, npoint, error = True)
r1B     = md.bin_up_data(r1, npoint)
r1erB   = md.bin_up_data(r1er, npoint, error = True)


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
	yy.append(r1B[i])

	xx.append(nhiB[i])
	yy.append(r1B[i])

	xx.append(nhiB[i]-dx2/2.)
	yy.append(r1B[i])

## Plot LH vs NHI ###
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

fig          = plt.figure(figsize=(16,10))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)

# major_xticks = np.arange(0., 500., 10.)
# minor_xticks = np.arange(0., 500., 10.)
major_yticks = np.arange(1., 11., 1.)                                              
minor_yticks = np.arange(1., 11., 0.5)


plt.errorbar(nh1, sig1, xerr=nh1*0., yerr=er1*0., color='b', marker='o', ls='None', markersize=mks, markeredgecolor='b', markeredgewidth=1, label='')
xerb1,  = plt.plot(nh1, sig1, color='b', marker='o', ls='None', markersize=mks, markeredgecolor='b', markeredgewidth=1, label='')

plt.errorbar(nh3, sig3, xerr=nh3*0., yerr=er3*0., color='purple', marker='o', ls='None', markersize=mks, markeredgecolor='purple', markeredgewidth=1, label='')
xerb2,  = plt.plot(nh3, sig3, color='purple', marker='o', ls='None', markersize=mks, markeredgecolor='purple', markeredgewidth=1, label='')

# plt.errorbar(nh2, sig2, xerr=nh2*0., yerr=er2, color='0', marker='o', ls='None', markersize=mks, markeredgecolor='k', markeredgewidth=1, label='')
xerb3,  = plt.plot(nh2, sig2, color='0', marker='o', ls='None', markersize=mks, markeredgecolor='k', markeredgewidth=1, label='')
plt.fill_between(nh2, sig2-er2, sig2+er2, color='lightgray', alpha=0.5)

plt.errorbar(nhi, r1, xerr=nhier, yerr=r1er, color='r', marker='d', ls='None', markersize=mks, markeredgecolor='r', markeredgewidth=1, label='')
xerb4,  = plt.plot(nhi, r1, color='r', marker='d', ls='None', markersize=mks, markeredgecolor='r', markeredgewidth=1, label='')

plt.plot([0.,500.], [r_av,r_av], 'k', mew=2, linewidth=2, linestyle='-.', marker='o', markerfacecolor='b', markersize=0, label='')
plt.plot([0.,500.], [pl_av,pl_av], 'k', mew=2, linewidth=2, linestyle='--', marker='o', markerfacecolor='b', markersize=0, label='')

# Bin-up #
# plt.plot(xx, yy, 'k', mew=2, linewidth=2, linestyle='-', marker='o', markerfacecolor='b', markersize=0, label='')
# plt.errorbar(nhiB, r1B, xerr=nhierB, yerr=r1erB, color='r', marker='o', ls='None', markersize=mks, markeredgecolor='r', markeredgewidth=1, label='')

plt.title('', fontsize=0)
plt.ylabel(r'$\mathrm{L_{H} = 4\pi \mathcal{R}/N_{H}\ [10^{-31}WH^{-1}]} $', fontsize=fts, fontweight='normal')
plt.xlabel(r'$\mathrm{N_{H} [10^{20} cm^{-2}]}$', fontsize=fts, fontweight='normal')

plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=15)
# plt.yscale('log')
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
plt.ylim(1.,10.)

axbox = ax.get_position()
leg   = plt.legend([xerb4, xerb1, xerb3, xerb2], [r'$\mathrm{This\ work}$', r'$\mathrm{X_{CO}}$$\mathrm{=1}$$\times$$\mathrm{10}$$^{20}$ $\mathrm{(PLC2014a)}$', \
	r'$\mathrm{X_{CO}}$$\mathrm{=2}$$\times$$\mathrm{10}$$^{20}$ $\mathrm{(PLC2014a)}$',\
	 r'$\mathrm{X_{CO}}$$\mathrm{=3}$$\times$$\mathrm{10}$$^{20}$ $\mathrm{(PLC2014a)}$'], \
	 fontsize=23, loc=(axbox.x0-0.11, axbox.y0+0.6), numpoints=1, handletextpad=-0.5)
leg.get_frame().set_linewidth(0.0)

plt.tight_layout()

plt.savefig('LH_vs_nh.eps', bbox_inches='tight', pad_inches=0.1, format='eps', dpi=600)
# plt.savefig('LH_vs_nh.png', bbox_inches='tight', pad_inches=0.03, format='png', dpi=100)
plt.show()
## END - PLOT ##