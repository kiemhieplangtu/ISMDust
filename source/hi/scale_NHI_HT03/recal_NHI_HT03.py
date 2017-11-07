import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust') # add folder of Class

from common.myImport import *
from scipy.io.idl    import readsav
from collections     import OrderedDict


##================= MAIN ========================##
map_file = const.DATPATH + 'oh/haslam408_dsds_Remazeilles2014.fits'
factor   = (408./1420.4)**2.8

## Haslam continuum
tb408    = hp.read_map(map_file,field=0, nest=False, hdu=1, h=False, verbose=False)
info79   = txtDat.read_79_info(fname='../data/79src_radec_lb.txt')
tbg21    = hpx.get_tb_haslam(tb408, info79)
inf408   = readsav('../../oh/data/tb_408.sav')

## Read fit params of HI components from Carl papers ##
info      = txtDat.read_HT03_params(fname = '../result/component_fit_params.txt', asarray=True)
src       = info['src']
frac      = info['frac']
cnm       = info['cnm']
nhi       = info['nhi']
err_tau   = info['err_tau']
err_tspin = info['err_tspin']
v0        = info['v0']
v0er      = info['err_v0']
tkmax     = info['tkmax']
err_frac  = info['err_frac']

del_v     = info['del_v']
err_del_v = info['err_del_v']
tau       = info['tau']
tspin     = info['tspin']
tb        = info['t_peak']
err_tb    = info['err_t_peak']


### Re-calc. N(HI) ###
with open('names.csv', 'w') as csvfile:
    fieldnames = ['Tpeak', 'er', 'tau', 'er', 'Vlsr', 'er', 'delV', 'er', 'Tspin', 'er', 'Tkmax', 'NHI', 'er', 'CNM', 'Frac', 'er', 'Name']
    writer     = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()

ret = {}
for i in range(len(src)):
	newNHI = 0.
	woc    = cnm[i]
	tauer  = err_tau[i]
	tser   = err_tspin[i]
	nhier  = 0.   
	if (woc==9):  # Warm component
		tser   = '----'
		woc    = '-'
		tauer  = '----'
		newNHI = const.HT03SCL*nhi[i]  ## already in 1e20
		wNHIer = const.HT03SCL*stats.get_NWNM_error(tb[i], err_tb[i], del_v[i], err_del_v[i])
		wNHIer = wNHIer*0.01  ## in Unit of 1e20
		string = '{:06.2f}   {:06.2f}   {:06.3f}   {:07.3f}   {:07.2f}   {:04.2f}   {:05.2f}   {:05.2f}   {:07.2f}   {:05.2f}   {:08.1f}   {:07.4f}   {:07.4f}   {}   {:06.2f}   {:06.2f}   {}'\
		.format(tb[i], err_tb[i], tau[i], err_tau[i], v0[i], v0er[i], del_v[i], err_del_v[i], tspin[i], err_tspin[i], tkmax[i], newNHI, round(wNHIer,4), cnm[i], frac[i], err_frac[i], src[i])
		print string
		xrow   = [tb[i], err_tb[i], tau[i], err_tau[i], v0[i], v0er[i], del_v[i], err_del_v[i], tspin[i], err_tspin[i], tkmax[i], newNHI, round(wNHIer,4), cnm[i], frac[i], err_frac[i], src[i]]
	else:
		etau   = np.exp(-tau[i])
		tbi    = tbg21[src[i]]

		newTs  = const.HT03SCL*tspin[i] + (const.HT03SCL-1.0)*tbi*etau/(1.0-etau)
		wTser  = stats.new_Tspin_error(tbi, const.HT03SCL, tser, tau[i], tauer)
		newNHI = const.NHICST*tau[i]*newTs*del_v[i]
		newNHI = 0.01*newNHI   ## in Unit of 1e20
		wNHIer = stats.get_NCNM_error(tau[i], tauer, newTs, wTser, del_v[i], err_del_v[i])
		wNHIer = wNHIer*0.01  ## in Unit of 1e20

		string = '{:06.2f}   {:06.2f}   {:06.3f}   {:07.3f}   {:07.2f}   {:04.2f}   {:05.2f}   {:05.2f}   {:07.2f}   {:05.2f}   {:08.1f}   {:07.4f}   {:07.4f}   {}   {:06.2f}   {:06.2f}   {}'\
		.format( tb[i], err_tb[i], tau[i], err_tau[i], v0[i], v0er[i], del_v[i], err_del_v[i], round(newTs,4), round(wTser,4), tkmax[i], round(newNHI,4), round(wNHIer,4), cnm[i], frac[i], err_frac[i], src[i] )
		xrow   = [ tb[i], err_tb[i], tau[i], err_tau[i], v0[i], v0er[i], del_v[i], err_del_v[i], round(newTs,4), round(wTser,4), tkmax[i], round(newNHI,4), round(wNHIer,4), cnm[i], frac[i], err_frac[i], src[i] ]
		print string

	row = OrderedDict()
	for i in range(len(fieldnames)):
		row[fieldnames[i]] = xrow[i]

	writer.writerow(row)

	if src[i] not in ret.keys():
		ret[src[i]] = {}

		ret[src[i]]['nhi']   = newNHI
		ret[src[i]]['nhier'] = wNHIer**2
	else:
		ret[src[i]]['nhi']   = ret[src[i]]['nhi']   + newNHI
		ret[src[i]]['nhier'] = ret[src[i]]['nhier'] + wNHIer**2


### Cal. Error - Squared-root ##
for sc in ret:
	ret[sc]['nhier'] = np.sqrt(ret[sc]['nhier'])

## print to: fit_params_recal.dat ##