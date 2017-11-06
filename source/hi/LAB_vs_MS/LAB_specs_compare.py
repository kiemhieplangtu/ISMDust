import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust') # add folder of Class

from common.myImport import *

## Calculate the uncertainties of factors #
 #
 # params 1D-array T
 #
 # return factor_uncertainties
 # 
 # Author Van Hiep ##
def cal_optthin_nhi(T):
	N  = len(T)
	dv = 1.03051546392
	s  = 0.
	for i in range(N):
		s = s + T[i]*dv

	vid1   = md.get_index(vLAB, 80.)
	vid2   = md.get_index(vLAB, 100.)
	stdTb  = np.std(tbLAB[vid1:vid2])

	thin   = s*1.824/100.
	thiner = md.uncertainty_of_WI(N, stdTb, dv)*1.824/100.

	return thin, thiner

###================= MAIN ========================###
## Infor for 78 MS sources
src79   = md.read_info_ms_79sc(fname = '../result/nhi_lb_79src_HT03.txt', asarray=True)
xl      = src79['l']
xb      = src79['b']
src     = src79['src']
thin    = src79['thin']
thiner  = src79['thiner']

msSpec  = md.read_hi_specs(fname = '../data/nhi_opac_specs.txt')

dirs    = os.listdir( '../data/LAB_specs/' )
delV    = 1.03051546392

for i in range(79):
	sc  = src[i]
	zl  = xl[i]
	zb  = xb[i]

	fle = sc+'.txt'
	if(fle not in dirs):
		continue

	vMS    = msSpec[sc]['v']
	texpMS = msSpec[sc]['Texp']

	vLAB,\
	tbLAB  = md.read_LAB_spec(fname = '../data/LAB_specs/'+sc+'.txt')
	

	# Plot spectra
	# fig    = plt.figure(figsize=(12,12))
	# ax     = fig.add_subplot(111); #ax.set_rasterized(True)

	# plt.plot(vMS, texpMS, 'k-', label='MS')
	# plt.plot(vLAB, tbLAB, 'r-', label='LAB')

	# plt.title(sc, fontsize = 35)
	# plt.ylabel('$T_{b} (K)$', fontsize = 35)
	# plt.xlabel('VLSR (km/s)', fontsize = 35)
	# plt.tick_params(axis='x', labelsize=20)
	# plt.tick_params(axis='y', labelsize=20)

	# plt.xlim(-100, 100)

	# plt.legend(loc='upper left', fontsize=18)

	# plt.savefig("figures/"+sc+'.eps', bbox_inches='tight', pad_inches=0.01, format='eps', dpi=60)
	# plt.close(fig)

	
thinLAB  = np.zeros(79)
thinLab  = np.zeros(79)
error    = np.zeros(79)
for i in range(79):
	sc  = src[i]
	zl  = xl[i]
	zb  = xb[i]

	fle = sc+'.txt'
	if(fle not in dirs):
		continue

	vMS    = msSpec[sc]['v']
	texpMS = msSpec[sc]['Texp']

	line   = open('../data/LAB_specs/'+sc+'.txt', "r").readlines()[2]
	nhiLAB = line.split(';')[4]
	nhiLAB = float(nhiLAB)/1e20


	vLAB,\
	tbLAB  = md.read_LAB_spec(fname = '../data/LAB_specs/'+sc+'.txt')
	nhiLab,\
	er     = cal_optthin_nhi(tbLAB)	

	thinLAB[i] = nhiLAB
	thinLab[i] = nhiLab
	error[i]   = er

# Plot - compare N*(HI)
fig    = plt.figure(figsize=(12,12))
ax     = fig.add_subplot(111); #ax.set_rasterized(True)

# plt.plot(thin, thinLab, 'rs', label='N*(HI) from Integral')
# plt.plot(thin, thinLAB, 'b.', label='N*(HI) from data file')

plt.errorbar(thin, thinLab, xerr=thiner, yerr=error, color='r', marker='s', ls='None', markersize=8, markeredgecolor='r', markeredgewidth=1, label='N*(HI) from Integral')
plt.errorbar(thin, thinLAB, xerr=thiner, yerr=error, color='b', marker='o', ls='None', markersize=2, markeredgecolor='b', markeredgewidth=1, label='N*(HI) from data file')

plt.plot([0,60], [0,60], 'k-')
plt.plot([0,60], [0,60*1.2], 'k:')

plt.title('LAB vs MS', fontsize = 35)
plt.ylabel('N(HI) - LAB', fontsize = 35)
plt.xlabel('N(HI) - MS', fontsize = 35)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)

plt.xlim(-1, 60)

plt.legend(loc='upper left', fontsize=18)
plt.show()


#############################
#         To fit            #
#############################
xdata     = thin
ydata     = thinLab
xerr      = thiner
yerr      = error

### Correlation:    ##
coxy      = md.cov_xy(xdata,ydata)
varx      = md.var(xdata)
vary      = md.var(ydata)
rho       = coxy/varx/vary

print ''
print '********* Pearson Coefficient *********'
print 'Pearson Coeff', md.pearson_coeff(xdata, ydata), ', ', rho
print 'Number of datapoints: ', len(xdata)
print ''

plt.plot(xdata/varx, ydata/vary, 'r*')
plt.xlabel('X/varx')
plt.ylabel('Y/vary')
plt.show()
### End - Correlation:##

########### ODR fit ############
xfit, yfit, mu, sig, m, ea = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=[1.1], xplot=[0., 60.], plot=False)
# xfit, yfit, mu, sig, m, ea, b, eb = md.do_linODRfit(xdata, ydata, xerr, yerr, lguess=[1.1, 0.1], xplot=[0., 60.], plot=False)

print 'Fit Results:'
print m
print ea

# print b
# print eb



## For plotting
mks    = 8
fts    = 36

plt.errorbar(thin, thinLab, xerr=thiner, yerr=error, color='k', marker='o', ls='None', markersize=8, markeredgecolor='k', markeredgewidth=1, capsize=0, label='$N^*_{HI}$l')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='None', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu-sig, mu+sig, color='0.5', alpha=0.5)

plt.plot([0,60], [0,60], 'k--')
# plt.plot([0,60], [0,60*1.2], 'k:')

plt.ylabel('$N^*_{HI}$ - LAB', fontsize=35)
plt.xlabel('$N^*_{HI}$ - HT03', fontsize=35)

plt.xlim(-1, 60)

plt.grid(True)
plt.tick_params(axis='x', labelsize=18)
plt.tick_params(axis='y', labelsize=15)
# plt.yscale('log')
# plt.xscale('log')
plt.legend(loc='upper left', fontsize=18)
plt.show()




# Plot for 
mpl.rcParams['axes.linewidth'] = 1.5
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rc('xtick', labelsize=15)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

fig          = plt.figure(figsize=(10,6))
ax           = fig.add_subplot(111); #ax.set_rasterized(True)                                 
major_xticks = np.arange(0., 60., 10.)
minor_xticks = np.arange(0., 60., 5.)
major_yticks = np.arange(0., 80., 10.)                                              
minor_yticks = np.arange(0., 80., 5.0)

plt.errorbar(thin, thinLab, xerr=thiner, yerr=error, color='k', marker='o', ls='None', markersize=6, markeredgecolor='k', markeredgewidth=1, capsize=0, label='')
plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='None', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu-sig, mu+sig, color='0.5', alpha=0.5)
plt.plot([0,60], [0,60], 'k--', label='$x = y$')

plt.xlabel(r'$\mathrm{N^*_{HI}\ (HT03)\ [10^{20} cm^{-2}]}$', fontsize=20)
plt.ylabel(r'$\mathrm{N^*_{HI}\ (LAB)\ [10^{20} cm^{-2}]}$', fontsize=20)

ax.set_xticks(major_xticks)                                                       
ax.set_xticks(minor_xticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
ax.tick_params(axis='x', labelsize=18, pad=5)
ax.tick_params(axis='y', labelsize=18)
ax.tick_params(which='both', width=1.5)
ax.tick_params(which='major', length=9)
ax.tick_params(which='minor', length=4)

if(0):
	ax.set_xlim(0.8, 50.)
	ax.set_ylim(0.6, 60.)

	ax.set_xscale('log')
	ax.set_yscale('log')

	plt.savefig('LAB_vs_HT03_logscale.eps', bbox_inches='tight', pad_inches=0.08, format='eps', dpi=600)

	plt.legend(loc='upper left', fontsize=18)

	plt.tight_layout()
	# plt.savefig('LAB_vs_HT03.eps', bbox_inches='tight', pad_inches=0.08, format='eps', dpi=600)

else:
	ax.set_xlim(-2.0, 50.)
	ax.set_ylim(-2.0, 60.)

	plt.legend(loc='upper left', fontsize=18)

	plt.tight_layout()
	# plt.savefig('LAB_vs_HT03.eps', bbox_inches='tight', pad_inches=0.08, format='eps', dpi=600)



# ax.set_xticks([5., 10., 20., 30.])
# ax.set_yticks([0.5, 1., 5., 10.])
# ax.get_xaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.0f'))
# ax.get_yaxis().set_major_formatter(mpl.ticker.FormatStrFormatter('%.1f'))

plt.show()



# Plot the ratios
plt.plot(xb, thinLab/thin, 'r.')
plt.show()