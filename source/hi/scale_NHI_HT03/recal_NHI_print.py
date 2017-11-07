import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust') # add folder of Class

from common.myImport import *

##================= MAIN ========================##

### 79 sources ###
ms79sc   = md.read_info_ms_79sc(fname = '../result/nhi_lb_79src_HT03.txt')
sc79     = ms79sc['src']
xl79     = ms79sc['l']
xb79     = ms79sc['b']
nhi03    = ms79sc['nhi']
nhier03  = ms79sc['nhier']
thin03   = ms79sc['thin']
thin03er = ms79sc['thiner']
cnm03    = ms79sc['cnm']
cnm03er  = ms79sc['cnmer']
wnm03    = ms79sc['wnm']
wnm03er  = ms79sc['wnmer']

old      = {}
for i in range(len(sc79)):
	if sc79[i] not in old.keys():
		old[sc79[i]]           = {}
		old[sc79[i]]['l']      = xl79[i]
		old[sc79[i]]['b']      = xb79[i]
		old[sc79[i]]['nhi']    = nhi03[i]
		old[sc79[i]]['nhier']  = nhier03[i]
		old[sc79[i]]['thin']   = thin03[i]
		old[sc79[i]]['thiner'] = thin03er[i]
		old[sc79[i]]['cnm']    = cnm03[i]
		old[sc79[i]]['cnmer']  = cnm03er[i]
		old[sc79[i]]['wnm']    = wnm03[i]
		old[sc79[i]]['wnmer']  = wnm03er[i]

## Read fit params of HI components from Carl papers ##
cols  = ['t_peak','err_t_peak','tau','err_tau','v0','err_v0','del_v','err_del_v', 'tspin', 'err_tspin', 'tkmax', 'nhi', 'nhier' , 'cnm', 'frac', 'err_frac', 'src']
fmt   = ['f',     'f',          'f',  'f',      'f', 'f',     'f',    'f',         'f',      'f',         'f',     'f',  'f',     'f',    'f',   'f',        's']
dat   = restore('fit_params_recal_78src.dat', 4, cols, fmt)
inf   = dat.read(asarray=True)

src       = inf['src']
frac      = inf['frac']
cnm       = inf['cnm']
nhi       = inf['nhi']
nhier     = inf['nhier']
err_tau   = inf['err_tau']
err_tspin = inf['err_tspin']
v0        = inf['v0']
v0er      = inf['err_v0']
tkmax     = inf['tkmax']
err_frac  = inf['err_frac']

del_v     = inf['del_v']
err_del_v = inf['err_del_v']
tau       = inf['tau']
tspin     = inf['tspin']
tb        = inf['t_peak']
err_tb    = inf['err_t_peak']

ret = {}
for i in range(len(src)):
	woc = cnm[i]
	if (woc==9):           # Warm component
		xnhi   = nhi[i]
		xnhier = nhier[i]
		xwnm   = nhi[i]
		xwnmer = nhier[i]
		xcnm   = 0.
		xcnmer = 0.
	else:                  # CNM
		xnhi   = nhi[i]
		xnhier = nhier[i]
		xwnm   = 0.
		xwnmer = 0.
		xcnm   = nhi[i]
		xcnmer = nhier[i]


	if src[i] not in ret.keys():
		ret[src[i]] = {}

		ret[src[i]]['nhi']   = xnhi
		ret[src[i]]['nhier'] = xnhier**2

		ret[src[i]]['cnm']   = xcnm
		ret[src[i]]['cnmer'] = xcnmer**2

		ret[src[i]]['wnm']   = xwnm
		ret[src[i]]['wnmer'] = xwnmer**2
	else:
		ret[src[i]]['nhi']   = ret[src[i]]['nhi']   + xnhi
		ret[src[i]]['nhier'] = ret[src[i]]['nhier'] + xnhier**2

		ret[src[i]]['wnm']   = ret[src[i]]['wnm']   + xwnm
		ret[src[i]]['wnmer'] = ret[src[i]]['wnmer'] + xwnmer**2

		ret[src[i]]['cnm']   = ret[src[i]]['cnm']   + xcnm
		ret[src[i]]['cnmer'] = ret[src[i]]['cnmer'] + xcnmer**2


### Cal. Error - Squared-root ##
for sc in ret:
	ret[sc]['nhier'] = np.sqrt(ret[sc]['nhier'])
	ret[sc]['cnmer'] = np.sqrt(ret[sc]['cnmer'])
	ret[sc]['wnmer'] = np.sqrt(ret[sc]['wnmer'])


### For 3C223 ###
ret['3C223']           = {}
ret['3C223']['nhi']    = 1.26*old['3C223']['nhi']
ret['3C223']['nhier']  = 1.26*old['3C223']['nhier']
ret['3C223']['thin']   = old['3C223']['thin']
ret['3C223']['thiner'] = old['3C223']['thiner']
ret['3C223']['cnm']    = 1.26*old['3C223']['cnm']
ret['3C223']['cnmer']  = 1.26*old['3C223']['cnmer']
ret['3C223']['wnm']    = 1.26*old['3C223']['wnmer']
ret['3C223']['wnmer']  = 1.26*old['3C223']['wnmer']

msnhi   = np.zeros(79)
newcl   = np.zeros(79)
msnhier = np.zeros(79)
newcler = np.zeros(79)
k       = 0
for sc in ret:
	string      = '{} \t {:10s} \t{:03.4f}   {:03.4f}\t{:03.4f}\t{:03.4f}\t{:03.4f}\t{:03.4f}\t{:03.4f}\t{:03.4f}\t{:03.4f}\t{:03.4f}'\
	.format(k, sc, old[sc]['l'], old[sc]['b'], ret[sc]['nhi'], ret[sc]['nhier'], 1.26*old[sc]['thin'], 1.26*old[sc]['thiner'], ret[sc]['cnm'], ret[sc]['cnmer'], ret[sc]['wnm'], ret[sc]['wnmer'] )
	print string                    ## print to: nhi_79src_recal.txt
	msnhi[k]    = old[sc]['nhi']
	newcl[k]    = ret[sc]['nhi']
	msnhier[k]  = old[sc]['nhier']
	newcler[k]  = ret[sc]['nhier']
	k += 1


print ret['3C223']

### Fit ###
print msnhi

xfit, yfit, mu, sig, m, ea, b, eb = md.do_linODRfit(msnhi, newcl, msnhier, newcler, lguess=[1.25, 0.5], plot=False)

### Plot ###
plt.figure(1, figsize=(16, 14))

plt.plot(xfit, mu, '-r', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

plt.errorbar(msnhi, newcl, xerr=msnhier, yerr=newcler, color='r', marker='o', ls='', markersize=6, markeredgecolor='r', markeredgewidth=1, label='data')

plt.plot(xfit, mu, '-k', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
plt.fill_between(xfit, mu-sig, mu+sig, color='0.5', alpha=0.5)

plt.plot([0,160],[0,160], 'k--', label='$x=y$')
plt.xlabel('N(HI) - HT03')
plt.ylabel('New - N(HI)')
plt.show()