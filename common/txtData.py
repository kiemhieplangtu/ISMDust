from basicImport import *
from restore     import restore

#######################################################################
# Block for reading data
# Author: Van Hiep
# Version: 08/2017
##################### - START - #######################################

## Read info of 30 SPONGE sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 12/2016
 # Author Van Hiep ##
def read_info_sponge_30src(fname = '../../oh/sponge/sub_data/30src_claire.txt', asarray=False):
	cols = ['src','l', 'b', 'cnm','cnm_er','wnm','wnm_er','nhi','nhi_er','thin','thin_er']
	fmt  = ['s',  'f', 'f', 'f',    'f',   'f',  'f',     'f',   'f',     'f',    'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=asarray)
	return dat

## Read info of 27 atomic sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_atomic_src(fname = '../../oh/result/27atomic_src.dat', asarray=False):
	cols = ['idx','src','l', 'b', 'nhi','nhier','thin','thiner', 'cnm','cnmer','wnm','wnmer']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=asarray)
	return dat

## Read info of 27 atomic sources with E(B-V) #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_ebv_atomic_src(fname = '../dust/ebv2nh/data/27atomic_src_with_ebv.dat', asarray=False):
	cols = ['idx','src','l', 'b', 'nhi','nhier','thin','thiner', 'cnm','cnmer','wnm','wnmer', 'ebv', 'ebver', 'av']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f',  'f',    'f',    'f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=asarray)
	return dat	

## Read NHI from 93src, all MS LOS recal. NHI #
 #
 # params string fname Filename
 # return dict info of N(HI)
 # 
 # version 5/2017
 # Author Van Hiep ##
def read_nhi_93src(fname = '../hi/result/nhi_thin_cnm_wnm_93src.txt', asarray=False):
	cols = ['indx', 'src', 'l', 'b', 'nhi', 'nhi_er', 'thin', 'thin_er', 'cnm', 'cnm_er', 'wnm', 'wnm_er']
	fmt  = ['i',     's',   'f','f',  'f',   'f',     'f',     'f',        'f',   'f',      'f',   'f'   ]
	data = restore(fname, 2, cols, fmt)
	return data.read(asarray=asarray)

## Read infor from 93src, all MS LOS recal. NHI, Thin, WCNM, CNM #
 #
 # params string fname Filename
 # return dict Info of NHI, Thin, WCNM, CNM
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_93src_info(fname = '../hi/result/nhi_thin_cnm_wnm_93src.txt'):
	from collections import OrderedDict

	cols   = ['indx', 'src', 'l', 'b', 'nhi', 'nhi_er', 'thin', 'thin_er', 'cnm', 'cnm_er', 'wnm', 'wnm_er']
	fmt    = ['i',     's',   'f','f',  'f',   'f',     'f',     'f',        'f',   'f',      'f',   'f'   ]
	data   = restore(fname, 2, cols, fmt)

	ms93sc = data.read()

	sc93   = ms93sc['src']
	id93   = ms93sc['indx']
	xl93   = ms93sc['l']
	xb93   = ms93sc['b']
	nhi    = ms93sc['nhi']
	nhier  = ms93sc['nhi_er']
	thin   = ms93sc['thin']
	thiner = ms93sc['thin_er']
	cnm    = ms93sc['cnm']
	cnmer  = ms93sc['cnm_er']
	wnm    = ms93sc['wnm']
	wnmer  = ms93sc['wnm_er']

	new    = OrderedDict()
	for i in range(len(sc93)):
		# if sc93[i] not in new.keys():
		new[sc93[i]]           = {}
		new[sc93[i]]['id']     = id93[i]
		new[sc93[i]]['l']      = xl93[i]
		new[sc93[i]]['b']      = xb93[i]
		new[sc93[i]]['nhi']    = nhi[i]
		new[sc93[i]]['nhier']  = nhier[i]
		new[sc93[i]]['thin']   = thin[i]
		new[sc93[i]]['thiner'] = thiner[i]
		new[sc93[i]]['cnm']    = cnm[i]
		new[sc93[i]]['cnmer']  = cnmer[i]
		new[sc93[i]]['wnm']    = wnm[i]
		new[sc93[i]]['wnmer']  = wnmer[i]

	return new

## Read NH and sigma353 from Planck #
 #
 # params string fname Filename
 # return dict info of NH and Sigma353
 # 
 # version 6/2017
 # Author Van Hiep ##
def read_planck_sigma_vs_nh(fname = '../dust/tau2nh/marc_data/sigma_e353_vs_N_H_Xco1.txt'):
	cols   = ['nh', 'sig', 'sd_sig']
	fmt    = ['f',  'f',    'f']
	data   = restore(fname, 0, cols, fmt)
	dat    = data.read()

	nh     = dat['nh']
	sig    = dat['sig']
	sd_sig = dat['sd_sig']
	return [nh, sig, sd_sig]

## Read info of E(B-V) for OH sources only #
 # src | l | b | E(B-V)SFD_online | err | E(B-V)S&F2011 | err
 #
 # params string fname Filename
 # return dict info
 # 
 # version 04/2017
 # Author Van Hiep ##
def read_ebv_for_oh_src(fname = 'data/ebv_sfd98_sf2011_for_oh_src.txt', sfd98=False):
	cols = ['idx','src','l', 'b', 'ebv','ebv_er', 'ebvsf','ebvsf_er', 'av']
	fmt  = ['i',  's',  'f', 'f', 'f',  'f'     , 'f',     'f'      , 'f' ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=True)

	if(sfd98):
		return dat['ebv'], dat['ebv_er'], dat['av'], dat['src']
	else:
		return dat['ebvsf'], dat['ebvsf_er'], dat['av'], dat['src']

## Read info of E(B-V) and Av#
 # src | l | b | E(B-V)SFD_online | err | E(B-V)S&F2011 | err | Av
 #
 # params string fname Filename
 # return dict info
 # 
 # version 08/2017
 # Author Van Hiep ##
def read_ebv_av(fname = 'data/ebv_sfd98_sf2011_for_93src.txt', sfd98=False):
	cols = ['src','l', 'b', 'ebv','ebv_er', 'ebvsf','ebvsf_er', 'av']
	fmt  = ['s',  'f', 'f', 'f',  'f'     , 'f',     'f'      , 'f' ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=True)

	if(sfd98):
		return dat['ebv'], dat['ebv_er'], dat['av'], dat['src']
	else:
		return dat['ebvsf'], dat['ebvsf_er'], dat['av'], dat['src']

## Read info of E(B-V) and Av for 93 src#
 # src | l | b | E(B-V)SFD_online | err | E(B-V)S&F2011 | err | Av
 #
 # params string fname Filename
 # return dict info
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_ebv_av_93src(fname = 'data/ebv_sfd98_sf2011_for_93src.txt', sfd98=False, asarray=True):
	cols = ['src','l', 'b', 'ebv','ebv_er', 'ebvsf','ebvsf_er', 'av']
	fmt  = ['s',  'f', 'f', 'f',  'f'     , 'f',     'f'      , 'f' ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray)

	return dat

## Read info of Av and E(B-V) for OH sources #
 # src | l | b | E(B-V)SFD_online | err | E(B-V)S&F2011 | err
 #
 # params string fname Filename
 # return dict info
 # 
 # version 07/2017
 # Author Van Hiep ##
def read_av_for_oh_src(fname = 'data/ebv_sfd98_sf2011_for_oh_src.txt'):
	cols  = ['idx','src','l', 'b', 'ebv','ebv_er', 'ebvsf','ebvsf_er', 'av']
	fmt   = ['i',  's',  'f', 'f', 'f',  'f'     , 'f',     'f'      , 'f' ]
	data  = restore(fname, 2, cols, fmt)
	dat   = data.read()

	ebv   = dat['ebvsf']
	ebver = dat['ebvsf_er']
	av    = dat['av']
	aver  = 3.1*np.array(ebver)
	src   = dat['src']

	ret   = {}
	for i in range(len(src)):
		if src[i] not in ret.keys():
			ret[src[i]] = {}

			ret[src[i]]['ebv']    = ebv[i]
			ret[src[i]]['ebv_er'] = ebver[i]
			ret[src[i]]['av']     = av[i]
			ret[src[i]]['aver']   = aver[i]

	return ret

## Read info of OH sources #
 # l,b, noh, noh_er
 #
 # params string fname Filename
 # return dict info
 # 
 # version 4/2017
 # Author Van Hiep ##
def read_noh(fname = '../../oh/result/total_noh65_21src.txt'):
	cols = ['idx','src','l', 'b', 'noh', 'noh_er']
	fmt  = ['i', 's',   'f', 'f', 'f',   'f']
	data = restore(fname, 2, cols, fmt)
	dat  = data.read()
	noh  = dat['noh']
	er   = dat['noh_er']
	src  = dat['src']
	l    = dat['l']
	b    = dat['b']

	ret  = {}
	for i in range(0,len(src)):
		# if dat['src'][i] not in ret.keys():
		ret[src[i]] = {}
		ret[src[i]]['noh']   = noh[i]
		ret[src[i]]['l']     = l[i]
		ret[src[i]]['b']     = b[i]
		ret[src[i]]['noher'] = er[i]

	return ret

## Read info of HI EM and ABS spectra #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 09/2017
 # Author Van Hiep ##
def read_hi_specs(fname = 'dark/hi/data/nhi_opac_specs.txt'):
	cols  = ['src','v', 'Texp', 'tau']
	fmt   = ['s',  'f', 'f',    'f'  ]
	data  = restore(fname, 3, cols, fmt)
	dat   = data.read(asarray=True)

	src   = dat['src']
	v     = dat['v']
	Texp  = dat['Texp']
	tau   = dat['tau']

	ret   = {}
	for i in range(len(src)):
		if src[i] not in ret.keys():
			ret[src[i]] = {}

			ret[src[i]]['v']    = [ v[i] ]
			ret[src[i]]['Texp'] = [ Texp[i]*0.5 ]
			ret[src[i]]['tau']  = [ tau[i] ]
		else:
			ret[src[i]]['v']    = ret[src[i]]['v']    + [ v[i] ]
			ret[src[i]]['Texp'] = ret[src[i]]['Texp'] + [ Texp[i]*0.5 ]
			ret[src[i]]['tau']  = ret[src[i]]['tau']  + [ tau[i] ]

	return ret


## Read info of 79 MS sources #
 # l,b, nhi, and nhi_error
 #
 # params string fname Filename
 # return dict info
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_ms_79src(fname = '../source/hi/data/79_src_nhi_iopscience.txt'):
	cols = ['src','h', 'mm', 's','d','min','ss', 'l','b','S','S_er', 'wnm', 'cnm', 'nhi']
	fmt  = ['s',  'f', 'f',  'f', 'f','f', 'f' , 'f','f','f',  'f', 'f',      'f',  'f' ]
	data = restore(fname, 32, cols, fmt)
	dat  = data.read()
	return dat

## Read LAB HI spectrum #
 #
 # params string fname Filename
 # return dict info
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_LAB_spec(fname = '../source/hi/data//LAB_specs/abs.txt'):
	vLAB, tbLAB, freq, wav = np.loadtxt(fname, skiprows=4, unpack=True)
	return vLAB, tbLAB

## Read info of 78 MS sources #
 # l,b, Ra, Dec
 #
 # params string fname Filename
 # return dict info
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_coord_ms_78src(fname = '../source/hi/data/78src_radec_lb.txt'):
	cols = ['l', 'b', 'src', 'ra', 'dec']
	fmt  = ['f', 'f',  's',  'f',  'f'  ]
	data = restore(fname, 3, cols, fmt)
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
def read_info_ms_79sc(fname = '../rearrange/nhi_lb_thin_78src.txt', asarray=False):
	cols = ['idx','src','l', 'b', 'nhi','nhier','thin','thiner', 'cnm','cnmer','wnm','wnmer']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=asarray)
	return dat

## Read infor from CSV #
 #
 # params string fname Filename
 # return dict info 
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_info_93src_csv(fname = '../doc/observed_src.csv',asarray=False):
	cols = ['idx', 'src', 'l', 'b', 'nhi','nhier','thin','thiner', 'cnm','cnmer','wnm','wnmer', 'ebv', 'ebver', 'av', 'ms', 'sp', 'hi', 'co', 'oh', 'noh', 'noher', 'tsig65', 'tsig67', 'mol', 'HIobs', 'semiAtom', 'cat1', 'cat2', 'cat3', 'ok', 'note']
	fmt  = ['i',    's',   'f', 'f', 'f',    'f',    'f',    'f',   'f',    'f',    'f',    'f', 'f',    'f',    'f', 'i',   'i',  'i',  'i',  'i',  'f', 'f',    'f',    'f',      'f',        'f',      'f',     'f',    'f',    'f',    'f',    's']
	data = restore(fname, 4, cols, fmt)
	dat  = data.readcsv(asarray=asarray)
	return dat

## Read infor from CSV file #
 #
 # params string fname Filename
 # params list cols Columns from data file
 # params list fmt  data-format of columns
 # params bool asarray  Read as array ?
 # return dict info 
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_info_csv(cols, fmt, fname = '../doc/observed_src.csv', skip=0, asarray=False):
	data = restore(fname, skip, cols, fmt)
	dat  = data.readcsv(asarray=asarray)
	return dat

## Get infor of 79 los #
 #
 # params str fname  File-name
 # params dict  info   Infor of the sources
 #
 # return list info
 # 
 # version 09/2017
 # Author Van Hiep ##
def read_79_info(fname='../data/79src_radec_lb.txt'):
	cols  = ['l','b','src','ra','dec']
	fmt   = ['f','f','s',  'f',  'f' ]
	dat   = restore(fname, 3, cols, fmt)
	info  = dat.read(asarray=True)

	return info

## Read fit params of HI components from Carl papers #
 #
 # params string fname Filename
 # params bool asarray  Read as array ?
 # return dict info 
 # 
 # version 10/2017
 # Author Van Hiep ##
def read_HT03_params(fname = '../result/component_fit_params.txt', asarray=False):
	## Read fit params of HI components from Carl papers ##
	cols   = ['t_peak','err_t_peak','tau','err_tau','v0','err_v0','del_v','err_del_v', 'tspin', 'err_tspin', 'tkmax', 'nhi', 'cnm', 'frac', 'err_frac', 'src']
	fmt    = ['f',     'f',          'f',  'f',      'f', 'f',     'f',    'f',         'f',      'f',         'f',     'f',   'i',    'f',   'f',        's']
	dat    = restore('../result/component_fit_params.txt', 3, cols, fmt)
	info   = dat.read(asarray=asarray)

	return info


#######################################################################
# End - Block for Reading data
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################