import sys, os
sys.path.insert(0, 'lib') # add folder of Class

import numpy             as np
from   restore           import restore


## Get the index of a given velocity #
 #
 # params list v-axis List of Velocity_axis
 # params float vel Velocity
 #
 # return int idx Index of vel in List of velocity_axis
 # 
 # Author Van Hiep ##
def get_vel_index(v_axis, vel):
    idx = (np.abs(np.array(v_axis)-vel)).argmin()
    return idx

## Get Vrange Indexes ##
 #
 # params 1-D-array v     VLSR
 # params float     lowv  Lower limit
 # params float     upv   Upper limit
 #
 # return list
 #
 # version 01/2017 
 # author Nguyen Van Hiep ##
def get_vrange_id(v, lowv, upv):
	vmax = get_vel_index(v, lowv)
	vmin = get_vel_index(v, upv)
	return [vmin, vmax]


#######################################################################
# Block for reading data
# Author: Van Hiep
# Version: 08/2017
##################### - START - #######################################



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
			ret[src[i]]['Texp'] = [ Texp[i] ]
			ret[src[i]]['tau']  = [ tau[i] ]
		else:
			ret[src[i]]['v']    = ret[src[i]]['v']    + [ v[i] ]
			ret[src[i]]['Texp'] = ret[src[i]]['Texp'] + [ Texp[i] ]
			ret[src[i]]['tau']  = ret[src[i]]['tau']  + [ tau[i] ]

	return ret


#######################################################################
# End - Block for Reading data
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################