from basicImport import *

#######################################################################
# Block for Correlation
# Author: Van Hiep
# Version: 08/2017
##################### - START - #######################################

## Cal. 1-sigma RMS of CO spec #
 #
 # params 1D-array v full-range VLSR
 # params 1D-array T full-range Tb
 # params float    v1, v2 vel-range to cal. RMS
 # params float    vmin, vmax good vel-range to cal. RMS
 #
 # return float rms
 #        1-D array v_fltr filtered vlsr 
 #        1-D array T_fltr filtered Tb 
 # 
 #  Version 6/2017
 # Author Van Hiep
 ##
def get_1sigma_spec(v,T, v1, v2, vmin=-50., vmax=50.):
	# fltr = np.extract([ (nh2>1.) & (nh2<4.) ], sig2)
	cond   = ( (v>vmin) & (v<v1) ) | ( (v>v2) & (v<vmax) )
	v_fltr = np.extract([ cond ], v)
	T_fltr = np.extract([ cond ], T)
	rms    = np.std(T_fltr)
	return rms, v_fltr, T_fltr

## Cal. Pearson correlation coefficient ##
 #
 # params list x x-data
 # params list y y-data
 #
 # return float r Pearson correlation coefficient
 #
 # version 04/2017 
 # author Nguyen Van Hiep ##
def pearson_coeff(x, y):
 	n   = len(x)
 	sxy = 0.
 	sx  = 0.
 	sy  = 0.
 	sx2 = 0.
 	sy2 = 0.
 	for i in range(n):
 		sxy = sxy + x[i]*y[i]
 		sx  = sx  + x[i]
 		sy  = sy  + y[i]
 		sx2 = sx2 + x[i]**2
 		sy2 = sy2 + y[i]**2

	t1 = n*sxy
	t2 = sx*sy

	t3 = n*sx2 - sx**2
	t3 = np.sqrt(t3)

	t4 = n*sy2 - sy**2
	t4 = np.sqrt(t4)

	r  = (t1-t2)/t3/t4

	return r

## Cal. covariance of 2 vectors ##
 #
 # params list x x-data
 # params list y y-data
 #
 # return float r Pearson correlation coefficient
 #
 # version 08/2017 
 # author Nguyen Van Hiep ##
def cov_xy(x, y):
 	n   = len(x)
 	mx  = np.mean(x)
 	my  = np.mean(y)

 	s   = 0.
 	for i in range(n):
 		s = s + (x[i]-mx)*(y[i]-my)

	return s/float(n-1)

## Cal. variance of 1 vector ##
 #
 # params list x x-data
 #
 # return float r Variance of a vectror
 #
 # version 08/2017 
 # author Nguyen Van Hiep ##
def var(x):
 	n   = len(x)
 	mx  = np.mean(x)

 	s   = 0.
 	for i in range(n):
 		s = s + (x[i]-mx)**2

	ret = s/float(n-1)
	return np.sqrt(ret)

## Calculate the uncertainties of a product a*b #
 #
 # params float a
 # params float b
 # params float aer
 # params float ber
 #
 # return float ret Uncertainty of a*b
 # 
 # Author Van Hiep ##
def uncertainty_of_product(a, b, aer, ber):
	d1 = aer*b
	d1 = d1*d1

	d2 = ber*a
	d2 = d2*d2

	d  = np.sqrt(d1+d2)

	return d


## Calculate the uncertainties of a product a*x+b #
 # 
 #
 # params float a
 # params float b
 # params float aer
 # params float ber
 #
 # return float ret Uncertainty of a*b
 # 
 # Author Van Hiep ##
def nh_uncert_from_proxies(x, xerr, a, aErr, b, bErr, corrEff=0.):
	d1 = x*aErr
	d1 = d1*d1

	d2 = bErr**2

	d3 = 2.0*x*corrEff*aErr*bErr

	d4 = a*xerr
	d4 = d4*d4

	r  = d1 + d2 + d3 + d4
	return np.sqrt(r)

## Calculate the uncertainties of ratio a/b #
 #
 # params float a
 # params float b
 # params float aer
 # params float ber
 #
 # return float ret Uncertainty of a/b
 # 
 # Author Van Hiep ##
def uncertainty_of_ratio(a, b, aer, ber):
	r  = np.abs(a/b)
	d1 = aer
	d1 = d1*d1

	d2 = r*ber
	d2 = d2*d2

	d  = np.abs(1.0/b) * np.sqrt(d1+d2)

	return d

## Calculate the uncertainties of a-b #
 #
 # params float ar
 # params float br
 #
 # return float er Uncertainty of a-b
 # 
 # Author Van Hiep ##
def uncertainty_of_diff(ar, br):
	e1 = ar*ar
	e2 = br*br
	er = np.sqrt(e1+e2)
	# print e1, e2, e1+e2-2.0*559.69e40
	return er

## Cal. uncertainty in the mean #
 #
 # params list lst List of numbers
 # return float ret Uncertainty in the mean
 # 
 # Author Van Hiep ##
def cal_uncertainty_in_mean(lst):
	n    = len(lst)
	mean = sum(lst)/float(n)

	s    = 0
	for i in range(n):
		s = s + (lst[i] - mean)**2

	s = np.sqrt(s)
	return s/n

## Calculate the uncertainties of factors #
 #
 # params 1D-array factr y-axis
 # params 1D-array nhi N(HI)
 # params 1D-array nhi_er Uncertainties of N(HI)
 # params 1D-array nh N(H)
 # params 1D-array nh_er Uncertainties of N(H)
 #
 # return factor_uncertainties
 # 
 # Author Van Hiep ##
def uncertainty_of_factors(factr, nhi, nhi_er, nh, nh_er):
	d1 = nhi_er/nhi
	d1 = d1**2

	d2 = nh_er/nh
	d2 = d2**2

	d  = np.sqrt(d1+d2)*factr

	return d

## Calculate the uncertainty of Integrated Intensity #
 #
 # params int/float N       Number of channels
 # params float     Sigma   Noise of Tb
 # params float     dv      Channel width
 #
 # return float     Uncertainty of Integrated Intensity
 # 
 # Author Van Hiep ##
def uncertainty_of_WI(N, sigma, dv):
	return np.sqrt(N)*sigma*dv


## Get  #
 #
 # params float tbg    background
 # params float A        scaling-factor
 # params float tser   Error of Tspin
 # params float tau      tau
 # params float tauer      Error of Tau
 #
 # return float Error of new Tspin
 # 
 # version 09/2017
 # Author Van Hiep ##
def new_Tspin_error(tbg, A, tser, tau, tauer):
 	d1 = A*tser
 	d1 = d1 * d1

 	zz = np.exp(-tau)
 	d2 = (A-1.0)*tbg*(-zz)*tauer/(1.0-zz)**2
 	d2 = d2 * d2

 	d  = np.sqrt(d1+d2)

 	return d

## Get  Error of CNM components#
 #
 # params float tau      tau
 # params float tauer    Error of Tau
 # params float Ts       Ts
 # params float Tser     Error of Tspin
 # params float wid      Width
 # params float wider    Error of Width
 #
 # return float Error of new N(HI)
 # 
 # version 09/2017
 # Author Van Hiep ##
def get_NCNM_error(tau, tauer, Ts, Tser, wid, wider):
	d1 = Ts*wid*tauer
	d1 = d1 * d1

	d2 = Ts*tau*wider
	d2 = d2 * d2

	d3 = tau*wid*Tser
	d3 = d3 * d3

	d  = np.sqrt(d1+d2+d3)

	return 1.93988*d #1e18

## Get  error of WNM component #
 #
 # params float Tb       Tb
 # params float Tber     Error of Tb
 # params float wid      Width
 # params float wider    Error of Width
 #
 # return float Error of new N(HI)  in unit of 1e18
 # 
 # version 09/2017
 # Author Van Hiep ##
def get_NWNM_error(tb, tber, wid, wider):
	d1 = tb*wider
	d1 = d1 * d1

	d2 = wid*tber
	d2 = d2 * d2

	d  = np.sqrt(d1+d2)

	return 1.93988*d  #1e18

#######################################################################
# End - Block for Correlation
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################