# Use python 3 consistent printing and unicode
from __future__ import print_function
from __future__ import unicode_literals

import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import healpy            as hp
import scipy             as sp
import inspect
import copy

from   restore           import restore
from   mpfit             import mpfit
from   scipy.odr         import *


#######################################################################
# Block for Block for functions for fitting
# Author: Van Hiep
# Version: 08/2017
##################### - START - #######################################

## Linear fucntion ##
 #
 # params list/float p Parameters
 # params 
 #
 # return 
 #
 # version 03/2017 
 # author Nguyen Van Hiep ##
def lin_fc(p, x):
	if(len(p) == 1):
		m    = p
		return m*x
	elif(len(p) == 2):
		m, b = p
		return m*x+b
	else:
		print("ODR - lin_fc - Error !")

## Linear fucntion ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # def tb_exp(x,y,bg,tau,v0,wid,tex,cont,err=1.):
 # author Nguyen Van Hiep ##
def lin_fcn(p, fjac=None, x=None, y=None, err=None):
	if (len(p) == 1):
		model = p[0]*x
	else:
		model = p[0] * x + p[1]

	status = 0
	return [status, (y-model) / err]

## Do linear fitting using MPFIT ##
 #
 # params list xdata X-data
 # params list ydata Y-data
 # params list xerr  X-errors
 # params list yerr  Y-errors
 # params list lguess Guesses
 #
 # return 1-D array xfit X-Fit
 #        1-D array yfit Y-Fit
 #        float mu  Mean value of Y-Fit
 #        float sig Standard-deviation of Y-Fit
 #        float m   slope
 #        float ea  Error of slope
 #
 # return options   
 #        float b   offset
 #        float eb  Error of offset
 #
 # version 08/2017 
 # author Nguyen Van Hiep ##
def do_linMPfit(xdata, ydata, xerr, yerr, lguess=[1.0], xplot=[], plot=True):
	if(len(xplot)==0):
		xplotMin = xdata.min()
		xplotMax = xdata.max()
	else:
		xplotMin = xplot[0]
		xplotMax = xplot[1]

 	npar    = len(lguess)
	guessp  = np.array(lguess, dtype='float64')
	plimd   = [[False,False]]*npar
	plims   = [[0.,0.]]*npar
	parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]}
	pname   = ['slope','offset']
	pfix    = [False]*npar

	parinfo = []
	for i in range(len(guessp)):
		parinfo.append(copy.deepcopy(parbase))

	for i in range(len(guessp)):
		parinfo[i]['value']   = guessp[i]
		parinfo[i]['fixed']   = pfix[i]
		parinfo[i]['parname'] = pname[i]
		parinfo[i]['limited'] = plimd[i]

	x  = xdata.astype(np.float64)
	y  = ydata.astype(np.float64)
	er = yerr.astype(np.float64)

	fa = {'x':x, 'y':y, 'err':er}
	mp = mpfit(lin_fcn, guessp, parinfo=parinfo, functkw=fa, quiet=True)

	# ## ********* Results ********* ##
	print ('********* Results from Linear MPFIT *********')
	abp   = mp.params
	abper = mp.perror
	for i in range(len(parinfo)):
		print(("  %s = {1:10.5g} +/- {2:10.5g}"+
           "          (Starting guess: {3:10.5g})").\
            format(parinfo[i]['parname'],abp[i],abper[i],lguess[i]))

	print('')

	if(npar == 1):
		a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
		xfit  = np.linspace(xplotMin, xplotMax, 20)
		yfit  = a[:, None] * xfit
		mu    = yfit.mean(0)
		sig   = 1.0*yfit.std(0)
		fit   = abp[0]*x

		m     = round(abp[0],4)
		ea    = round(abper[0],4)

		return xfit, yfit, mu, sig, m, ea

	else:

		a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
		b     = np.array([ abp[1]-abper[1], abp[1]+abper[1] ])
		xfit  = np.linspace(xplotMin, xplotMax, 20)
		yfit  = a[:, None] * xfit + b[:, None]
		mu    = yfit.mean(0)
		sig   = 1.0*yfit.std(0)
		fit   = abp[0]*x+abp[1]

		m     = round(abp[0],4)
		b     = round(abp[1],4)
		ea    = round(abper[0],4)
		eb    = round(abper[1],4)

		return xfit, yfit, mu, sig, m, ea, b, eb

	########### Plot Example ############
	# plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
	# plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='MPFIT linear fit')
	# plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

	# plt.title('', fontsize=30)
	# plt.xlabel('$Radiance\ . 10^{11} (Wcm^{-2}sr^{-1})$', fontsize=35)
	# plt.ylabel('$N_{HI}/10^{20} (cm^{-2}$)', fontsize=35)
	# # plt.xlim(0.0, 2.0)
	# # plt.ylim(-1.0, 6.0)
	# plt.grid(True)
	# plt.tick_params(axis='x', labelsize=18)
	# plt.tick_params(axis='y', labelsize=15)

	# # plt.text(0.001, 22., '$Fit: N_{H} = ['+str(m)+'\pm'+str(ea) +']\cdot10^{31} R$', color='blue', fontsize=20)
	# plt.text(0.02, 22., '$fit = ['+str(m)+'\pm'+str(ea) +']\cdot10^{31} R + ['+str(b)+'\pm'+str(eb)+']\cdot10^{31}$', color='blue', fontsize=20)
	# plt.legend(loc='upper left', fontsize=18)
	# # plt.savefig("test.png",bbox_inches='tight')
	# # for i in range(len(src)):
	# # 	if(oh[i] > 0):
	# # 		plt.annotate('('+str(src[i])+', '+ str(xl[i])+', '+str(xb[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
	# # 	            xytext=(-50.,30.), textcoords='offset points',
	# # 	            arrowprops=dict(arrowstyle="->"),fontsize=12,
	# # 	            )
	# plt.show()
	########### END - MPFIT ############


## Do linear fitting using ODR ##
 #
 # params list xdata X-data
 # params list ydata Y-data
 # params list xerr  X-errors
 # params list yerr  Y-errors
 # params list lguess Guesses
 # params bool plot Plot or not
 #
 # return 1-D array xfit X-Fit
 #        1-D array yfit Y-Fit
 #        float mu  Mean value of Y-Fit
 #        float sig Standard-deviation of Y-Fit
 #        float m   slope
 #        float ea  Error of slope
 #
 # return options   
 #        float b   offset
 #        float eb  Error of offset
 #
 # version 08/2017 
 # author Nguyen Van Hiep ##
def do_linODRfit(xdata, ydata, xerr, yerr, lguess=[1.0], xplot=[], plot=True):
	x_uncertainties_ignored = False

	if(len(xplot)==0):
		xplotMin = xdata.min()
		xplotMax = xdata.max()
	else:
		xplotMin = xplot[0]
		xplotMax = xplot[1]

 	npar    = len(lguess)
	guessp  = np.array(lguess, dtype='float64')
	plimd   = [[False,False]]*npar
	plims   = [[0.,0.]]*npar
	parbase = {'value': 0., 'fixed': 0, 'parname': '', 'limited': [0, 0], 'limits': [0., 0.]}
	pname   = ['slope','offset']
	pfix    = [False]*npar

	parinfo = []
	for i in range(len(guessp)):
		parinfo.append(copy.deepcopy(parbase))

	for i in range(len(guessp)):
		parinfo[i]['value']   = guessp[i]
		parinfo[i]['fixed']   = pfix[i]
		parinfo[i]['parname'] = pname[i]
		parinfo[i]['limited'] = plimd[i]

	x  = xdata.astype(np.float64)
	y  = ydata.astype(np.float64)
	er = yerr.astype(np.float64)

	# Create a model for Orthogonal distance regression (ODR) fitting.
	lin_model = Model(lin_fc)
	# Create a RealData object using our initiated data from above.
	data      = RealData(x=xdata, y=ydata, sx=xerr, sy=yerr)
	# Set up ODR with the model and data.
	odr       = ODR(data, lin_model, beta0=lguess)
	# Run the regression.
	out       = odr.run()

	## ********* Results ********* ##
	abp       = out.beta
	abper     = out.sd_beta

	cov       = out.cov_beta   # parameter covariance matrix
	#"Quasi-chi-squared" is defined to be the [total weighted sum of squares]/dof
	#	i.e. same as numpy.sum((residual/adjusted_err)**2)/dof or
	#       numpy.sum(((out.xplus-x)/xerr)**2
	#                                +((y_data-out.y)/yerr)**2)/dof
	#	Converges to the conventional chi-squared for zero x uncertainties.
	quasi_chisq = out.res_var

	if quasi_chisq < 1.0 :
	    abper = abper/np.sqrt(quasi_chisq)
	if False: # debugging print statements
	    print("sd_beta",output.sd_beta)
	    print("cov_beta",output.cov_beta)
	    print("delta",output.delta)
	    print("eps",output.eps)
	    print("res_var",output.res_var)
	    print("rel_error",output.rel_error)

    ############################## PRINT THE RESULTS #############################
	print("***********************************************************")
	print("               ORTHOGONAL DISTANCE REGRESSION")
	print("***********************************************************\n")
	print("ODR algorithm stop reason: " + out.stopreason[0])
	print("\nFit {0} Data points: ".format(len(xdata)))
	print("To Model :")
	print(" ",inspect.getsource(lin_fc))
	if x_uncertainties_ignored:
	    print("** WARNING: x uncertainties set to zero in fit. **\n")

	print("Estimated parameters and uncertainties")
	for i in range(len(abp)) :
	    print(("  %s = {1:10.6g} +/- {2:10.6g}"+
           "          (Starting guess: {3:10.6g})").\
            format(parinfo[i]['parname'],abp[i],abper[i],lguess[i]))

	# number of degrees of freedom for fit
	dof = len(xdata) - len(abp)

	# Note that odr, unlike curve_fit, returns the standard covariance matrix.
	print("\nStandard Covariance Matrix : \n", cov, "\n")

	print("\nCorrelation Matrix :")
	CorrCoeff = 1.
	for i,row in enumerate(cov):
	    for j in range(len(abp)) :
	        if (i != j):
	            CorrCoeff = cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])
	        print("{0:< 8.3g}".format(cov[i,j]/np.sqrt(cov[i,i]*cov[j,j])),
	                end="")
	            # Newbie Notes: "{0:< 8.3g}" left justifies output with space in
	            #   front of positive numbers, with 3 sig figs;
	            #   end="" suppresses new line.
	    print()
	print("\nCorrelation Coeff. of Parameters (m,b) : \n", CorrCoeff, "\n")
	# Calculate initial residuals and the 'adjusted error' for each data point
	delta = out.delta  # estimated x-component of the residuals
	eps   = out.eps    # estimated y-component of the residuals
	# (xstar,ystar) is the point where the 'residual line' (in black)
	#   intersects the 'ellipse' created by xerr & yerr.
	xstar        = ( xerr*np.sqrt( ((yerr*delta)**2) /
	                ( (yerr*delta)**2 + (xerr*eps)**2 ) ) )
	ystar        = ( yerr*np.sqrt( ((xerr*eps)**2) /
	                ( (yerr*delta)**2 + (xerr*eps)**2 ) ) )
	adjusted_err = np.sqrt(xstar**2 + ystar**2)
	# residual is positive if the point lies above the fitted curve,
	#             negative if below
	residual     = np.sign(ydata-lin_fc(abp,xdata))*np.sqrt(delta**2 + eps**2)

	print("\nQuasi Chi-Squared/dof   = {0:10.5f}".format(quasi_chisq) )
	print("   WARNING:Above CDF is not valid for large x uncertainties!")

	print("\nTo use Monte Carlo simulation to more accurately estimate CDF for")
	print('      large x uncertainties, re-run program with ')
	print('     "Run_Monte_Carlo_CDF_Estimator = True" and')
	print('     "Number_of_MC_iterations >= 1000." This may take some time\n')

	if(plot):
		# create figure with light gray background
		fig = plt.figure(facecolor="0.98")

		# 3 rows, 1 column, subplot 1
		#   3 rows are declared, but there are only 2 plots; this leaves room for text
		#       in the empty 3rd row
		fit = fig.add_subplot(211)
		# remove tick labels from upper plot (for clean look)
		fit.set_xticklabels( () )
		plt.ylabel("Y")

		plt.title("Orthogonal Distance Regression Fit to Data")
		# Plot data as red circles, and fitted function as (default) line.
		#   For a smooth look,generate many x values for plotting the model
		x_model = np.arange( min(xdata),max(xdata),
		                        (max(xdata)-min(xdata))/1000.)
		fit.plot(xdata,ydata,'ro', x_model, lin_fc(abp,x_model))
		# Add error bars on data as red crosses.
		fit.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, fmt='r+')
		fit.set_yscale('linear')
		#   draw starting guess as dashed green line ('r-')
		fit.plot(x_model, lin_fc(lguess,x_model), 'g-', label="Start", linestyle="--")

		a = np.array([out.xplus,xdata])   # output.xplus = x + delta
		b = np.array([out.y,ydata])       # output.y = f(p, xfit), or y + epsilon
		fit.plot(np.array([a[0][0],a[1][0]]), np.array([b[0][0],b[1][0]]),
		         'k-', label = 'Residuals')
		for i in range(1,len(ydata)):
		    fit.plot(np.array([a[0][i],a[1][i]]), np.array([b[0][i],b[1][i]]),
		         'k-')
		fit.legend(loc='upper left')
		fit.grid()

		# separate plot to show residuals
		residuals = fig.add_subplot(212) # 3 rows, 1 column, subplot 2
		residuals.errorbar(x=xdata,y=residual,yerr=adjusted_err,
		                   			fmt="r+", label = "Residuals")
		# make sure residual plot has same x axis as fit plot
		residuals.set_xlim(fit.get_xlim())
		# Draw a horizontal line at zero on residuals plot
		plt.axhline(y=0, color='b')
		# Label axes
		plt.xlabel("X")
		# These data look better if 'plain', not scientific, notation is used, and if
		#   the tick labels are not offset by a constant (as is done by default).
		plt.ticklabel_format(style='plain', useOffset=False, axis='x')
		plt.ylabel("Residuals")

		residuals.grid()

		plt.show()

	##### END ####
	if(npar == 1):
		a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
		xfit  = np.linspace(xplotMin, xplotMax, 20)
		yfit  = a[:, None] * xfit
		mu    = yfit.mean(0)
		sig   = 1.0*yfit.std(0)
		fit   = abp[0]*x

		m     = round(abp[0],4)
		ea    = round(abper[0],4)

		return xfit, yfit, mu, sig, m, ea

	else:

		a     = np.array([ abp[0]-abper[0], abp[0]+abper[0] ])
		b     = np.array([ abp[1]-abper[1], abp[1]+abper[1] ])
		xfit  = np.linspace(xplotMin, xplotMax, 20)
		yfit  = a[:, None] * xfit + b[:, None]
		mu    = yfit.mean(0)
		sig   = 1.0*yfit.std(0)
		fit   = abp[0]*x+abp[1]

		m     = round(abp[0],4)
		b     = round(abp[1],4)
		ea    = round(abper[0],4)
		eb    = round(abper[1],4)

		return xfit, yfit, mu, sig, m, ea, b, eb


	########### Plot Example ############
	# plt.errorbar(xdata, ydata, xerr=xerr, yerr=yerr, color='r', marker='o', ls='None', markersize=8, markeredgecolor='b', markeredgewidth=1, label='data')
	# plt.plot(xfit, mu, '-b', mew=2, linewidth=2, linestyle='solid', marker='o', markerfacecolor='b', markersize=0, label='ODR linear fit')
	# plt.fill_between(xfit, mu - sig, mu + sig, color='0.5', alpha=0.5)

	# plt.title('', fontsize=30)
	# plt.xlabel('$Radiance\ [10^{-11} Wcm^{-2}sr^{-1}]$', fontsize=35)
	# plt.ylabel('$N_{HI} [10^{20} cm^{-2}]$', fontsize=35)
	# plt.xlim(0.2, 10.0)
	# # plt.ylim(-1.0, 6.0)
	# plt.grid(True)
	# plt.tick_params(axis='x', labelsize=18)
	# plt.tick_params(axis='y', labelsize=15)
	# # plt.ticklabel_format(axis='x', style='sci', scilimits=(0,0))
	# plt.yscale('log')
	# plt.xscale('log')

	# plt.text(0.21, 12., '$N_{H} = ['+str(m)+'\pm'+str(ea) +']\cdot10^{31} R + ['+str(b)+'\pm'+str(eb)+']\cdot10^{20}$', color='blue', fontsize=20)
	# plt.text(0.21, 10, '(Data points with no CO and OH detected)', color='k', fontsize=20)
	# # plt.text(0.001, 12., '$Fit: N_{H} = ['+str(m)+'\pm'+str(ea) +']\cdot10^{31} R$', color='blue', fontsize=20)
	# # plt.gca().set_xscale('log', basex=10)
	# # plt.gca().set_yscale('log', basex=10)
	# plt.legend(loc='upper left', fontsize=18)
	# # for i in range(len(src)):
	# # 	if(oh[i] > 0):
	# # 		plt.annotate('('+str(src[i])+', '+ str(xl[i])+', '+str(xb[i])+')', xy=(xdata[i], ydata[i]), xycoords='data',
	# # 	            xytext=(-50.,30.), textcoords='offset points',
	# # 	            arrowprops=dict(arrowstyle="->"),fontsize=12,
	# # 	            )
	# plt.show()
	########### END - ODR ############


#######################################################################
# End - Block for functions for fitting
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################










#######################################################################
# End - Block for List
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################

## Find common elements of 2 lists ##
 #
 # params list x 
 # params list y 
 #
 # return list List of common elements
 #
 # version 12/2017 
 # author Nguyen Van Hiep ##
def common_elements(x, y):
 	return list(set(x) & set(y))

## Find Different elements of 2 lists ##
 #
 # params list x 
 # params list y 
 #
 # return list List of diff elements
 #
 # version 12/2017 
 # author Nguyen Van Hiep ##
def diff_elements(x, y):
 	return list(set(x) - set(y))

#######################################################################
# End - Block for List
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################











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

#######################################################################
# End - Block for Correlation
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################


#######################################################################
# Block for Placnk whole-map data
# Author: Van Hiep
# Version: 08/2017
##################### - START - #######################################


## Cal E(B-V) from Planck 2014 R1.2 #
 #
 # params array tau_map  Map of Tau353
 # params array tauer_map  Map of Tau353_error
 # params array ebv_map  Map of E(B-V)
 # params dict  info     Infor of the sources
 #
 # return list info
 # 
 # version 08/2017
 # Author Van Hiep ##
def cal_ebv_from_planckR12(tau_map, tauer_map, ebv_map, info):
	# Define constants #
	deg2rad = np.pi/180.

	fct     = 1.49e4
	fct_er  = 0.03e4

	## sources
	src    = info['src']  ## src
	xl     = info['l']
	xb     = info['b']
	n      = len(src)

	nside  = hp.get_nside(ebv_map)
	res    = hp.nside2resol(nside, arcmin=False)

	# OK - Go #
	ebv    = np.zeros(n)
	ebver  = np.zeros(n)
	for i in range(n):
		theta = (90.0-xb[i])*deg2rad
		phi   = xl[i]*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (ebv_map[pix] > -0.000001) : # Some pixels not defined
			xtau   = tau_map[pix]
			xtauer = tauer_map[pix]
			val    = ebv_map[pix]

			d1     = fct_er/fct
			d2     = xtauer/xtau
			dr     = np.sqrt(d1**2 + d2**2)
			err    = val*dr

		ebv[i]   = val
		ebver[i] = err

	return ebv, ebver


## Get tau353 and its error #
 #
 # params array tauMap  Map of tau353
 # params array errMap  Error map of tau353
 # params dict  info    Infor of the sources
 #
 # return 1-D arrays tau353 and tauer
 # 
 # version 12/2016
 # Author Van Hiep ##
def get_tau(tauMap, errMap, info):
	# src and l,b
	src = info['src']
	xl  = info['l']
	xb  = info['b']
	n   = len(src)

	# Define constants
	deg2rad = np.pi/180.
	nside   = hp.get_nside(tauMap)
	res     = hp.nside2resol(nside, arcmin=False)

	# OK - Go #
	tau353 = np.zeros(n)
	tauer  = np.zeros(n)
	for i in range(n):
		# Find the values of Tau353 and Err_tau353
		theta = (90.0-xb[i])*deg2rad
		phi   = xl[i]*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if ( (errMap[pix] >= 6.9e-11) and (errMap[pix] <= 0.00081) and (tauMap[pix] > -1.0e30) ): # Checked Invalid error & Some pixels not defined
			tau353[i] = tauMap[pix] 
			tauer[i]  = errMap[pix]
		else:
			print("Error! Error!")
			sys.exit()

	return tau353, tauer


## Get Radiance #
 #
 # params 1Darray tauMap     Map of tau353
 # params 1Darray tauErrMap  Map of tau353_error
 # params 1Darray rMap       Map of Radiance
 # params dict    info       Infor of the sources
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def get_radiance(tauMap, tauErrMap, rMap, info):
	# src and l,b
	src = info['src']
	xl  = info['l']
	xb  = info['b']
	n   = len(src)

	# Define constants #
	deg2rad = np.pi/180.

	# Constants to cal. error of Radiance from Tau353 (see PLC2014a for details)
	fct     = 0.0276
	fct_er  = 0.00072

	nside   = hp.get_nside(rMap)
	res     = hp.nside2resol(nside, arcmin=False)

	# OK - Go #
	r       = np.zeros(n)
	rer     = np.zeros(n)
	for i in range(n):
		theta = (90.0-xb[i])*deg2rad
		phi   = xl[i]*deg2rad
		pix   = hp.ang2pix(nside, theta, phi, nest=False)

		if (rMap[pix] > -0.000001) : # Some pixels not defined
			xtau   = tauMap[pix]
			xtauer = tauErrMap[pix]
			val    = rMap[pix]

			d1     = fct_er/fct
			d2     = xtauer/xtau
			dr     = np.sqrt(d1**2 + d2**2)
			err    = val*dr           ## Error of Radiance

		r[i]   = val
		rer[i] = err

	return r, rer


## Bin-up the data #
 #
 # params 1Darray  x  (In ascending order)    
 # params 1Darray  xer   
 # params 1Darray  y      
 # params 1Darray  yer   
 # params Int      n  Number of datapoints in a bin   
 #
 # return lists x, xer, y, yer
 # 
 # version 10/2017
 # Author Van Hiep ##
def bin_up_data(x, n, error=False):
	N     = len(x)
	mod   = N%n
	floor = N//n

	chnl  = floor
	if(mod != 0):
		chnl = chnl + 1

	indx = 0
	lst  = range(0,N,n)
	lmax = max(lst)

	xx   = np.zeros(chnl)
	for i in lst:
		nbin = n
		if( (i == lmax) and (mod != 0)) :
			nbin = mod

		if(error):
			xtemp = np.square(x[i:(i+nbin)])
			xtemp = np.sum(xtemp)
			xtemp = np.sqrt(xtemp)/nbin
		else:
			xtemp = np.mean(x[i:(i+nbin)])

		xx[indx] = xtemp
		indx     = indx + 1

	return xx


## Cal N(H) from tau353 #
 #
 # params array tau_map  Map of tau353
 # params array err_map  Error map of tau353
 # params dict  info     Infor of the sources
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def get_nh_from_tau353(tauMap, errMap, info):
	# Define constants #
	deg2rad     = np.pi/180.
	fukui_cf    = 2.10e26
	fk_fact_err = 0.0 #unknown

	## find Planck Conversion Factor (Dust opacity and Its error) ## 
	a, aer  = [1.39249e6, 4.9102e4] #6.6e-27, 0.66e-26, lowNHI
	b, ber  = [0., 0.]              ## If having intercerpt

	# Cal. tau353
	t353, t353er = get_tau(tauMap, errMap, info)

	# Calculate the NH from Planck factor #
	nh   = t353*a+b
	nher = nh_uncert_from_proxies(t353, t353er, a, aer, b, ber, corrEff=0.)

	return nh, nher, t353, t353er

## Cal N(H) from E(B-V) #
 #
 # params 1D array ebv     E(B-V)
 # params 1D array ebver   Error
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def get_nh_from_ebv(ebv, ebver):
	# a,aErr = [5.8e21, 0.0]        ## From Bohlin 1978   
	a, aer  = [113.912, 4.18695]    ## From S&F2011
	b, ber  = [0., 0.]

	nh      = ebv*a+b
	nher    = nh_uncert_from_proxies(ebv, ebver, a, aer, b, ber, corrEff=0.)

	return nh, nher, ebv, ebver

## Cal NH from Radiance #
 #
 # params array r_map    Map of Radiance
 # params dict  info     Infor of the sources
 #
 # return list info
 # 
 # version 04/2017
 # Author Van Hiep ##
def get_nh_from_radiance(tauMap, tauErrMap, rMap, info):
	a, aer = [4.65435e11, 0.170979e11]   ## N(H) = a.e31.R + b, NH/e20 = a.e11.R + b
	b, ber = [0., 0.]

	R, Rer = get_radiance(tauMap, tauErrMap, rMap, info)
	R      = R*1e-4     ## radiance W/m2/sr -> w/cm2/sr, xdata*1e11 => 1e7
	Rer    = Rer*1e-4

	nh     = R*a+b
	nher   = nh_uncert_from_proxies(R, Rer, a, aer, b, ber, corrEff=0.)

	return nh, nher, R, Rer

#######################################################################
# End - Block for Placnk whole-map data
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################



#######################################################################
# Block for Writing to file
# Author: Van Hiep
# Version: 08/2017
##################### - START - #######################################

## Write to file ##
 #
 # params str filename
 # params str string to write to file
 #
 # Example of a string ##
 # string = str(n)+ '\t'+ sc+ '\t'+ str( round(noh[sc]['l'],8) )+ '\t' + '\n'
 #
 # return Void
 #
 # version 08/2017 
 # author Nguyen Van Hiep ##
def write2file(filename, string):
	file = open(filename,'w') 
	file.write( string )
	file.close()

## Write to screen ##
 #
 # params str  style  Format to write
 # params list cols   Columns to write
 #
 # Example of a string ##
 # print '{:10s}{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}'\
 # 	.format(src, err, nh_er, nh2_er, nhier, n_h, nhi, nh2, nh_er, nh_er, nhier, nh2_er)
 #
 # return Void
 #
 # version 08/2017 
 # author Nguyen Van Hiep ##
def printCols(style='{:10s}{:08.4f}\t{:08.4f}', cols=[]):
	print (style)
	print (cols)
	# print '\''+style+'\''.format(cols)

#######################################################################
# End - Block for Writing to file
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################



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
	cols = ['idx','src','l', 'b', 'nhi','nhi_er','thin','thin_er', 'cnm','cnm_er','wnm','wnm_er']
	fmt  = ['i',  's',  'f', 'f', 'f',   'f',     'f',    'f'    , 'f',   'f',     'f',  'f'    ]
	data = restore(fname, 2, cols, fmt)
	dat  = data.read(asarray=asarray)
	return dat

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


#######################################################################
# End - Block for Reading data
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################


## Get x,y ID ranges for WCO contours #
 #
 # params float x1 X-min
 # params float x2 X-max
 # params float y1 Y-min
 # params float y2 Y-max
 #
 # return list of ID x-y ranges
 # 
 # version 9/2017
 # Author Van Hiep ##
def get_xy_id_ranges(x1=-180.,x2=180.,y1=-80.,y2=89.):
	x = np.arange(180.,-180.125,-0.125)
	y = np.arange(-80.,89.125,0.125)

	# xid1, = np.where(x==x2)
	# xid2, = np.where(x==x1)
	# yid1, = np.where(y==y1)
	# yid2, = np.where(y==y2)

	xid1 = (np.abs(x-x2)).argmin()
	xid2 = (np.abs(x-x1)).argmin()
	yid1 = (np.abs(y-y1)).argmin()
	yid2 = (np.abs(y-y2)).argmin()

	print('x1,x2,y1,y2', x1,x2,y1,y2)

	print('xid1, xid2, yid1, yid2', xid1, xid2, yid1, yid2)

	return yid1, yid2, xid1, xid2

## WCO map from Dame #
 #
 # params 
 # return 2-D array X, Y, Z
 # 
 # version 9/2017
 # Author Van Hiep ##
def Wco_map():
	from   astropy.io import fits

	image_file = os.getenv("HOME") + '/hdata/co/mlat_Wco_mom.fits'
	hdulist    = fits.open(image_file)

	# print hdulist[0].header
	data       = hdulist[0].data[:, ::-1]
	# data     = hdulist[0].data
	print('WCO data shape', data.shape )

	aa   = data[:,960:2881]
	bb   = data[:, 0:960]
	cc   = np.concatenate((aa, bb), axis=1)
	data = cc	


	Nx   = 2881
	Ny   = 1353

	x    = np.arange(180.,-180.125,-0.125)
	y    = np.arange(-80.,89.125,0.125)
	X, Y = np.meshgrid(x, y)

	return X,Y,data

## Patch of WCO map #
 #
 # params 2D-Arr X X-array
 # params 2D-Arr Y X-array
 # params 2D-Arr Z Z-array
 # params float x1 x-min
 # params float x2 x-max
 # params float y1 y-min
 # params float y2 y-max
 #
 # return list of ID x-y ranges
 # 
 # version 9/2017
 # Author Van Hiep ##
def patch_Wco_map(X,Y,Z, xmin=-180., xmax=180., ymin=-80., ymax=89.):
	# X,Y,Z       = Wco_map()
	x1,x2,y1,y2 = get_xy_id_ranges(xmin, xmax, ymin, ymax)

	xd = X[x1:x2,y1:y2]
	yd = Y[x1:x2,y1:y2]
	zd = Z[x1:x2,y1:y2]

	return xd, yd, zd


## Retreive a SINGLE value of 408 t_b from haslam et al. ##
 #
 # params float ell Galactic-longitude
 # params float bee Galactic-latitude
 #
 # return float Tbg_408 Background-Temperature at (l,b)
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_tb_408(ell,bee,tb_408):
	iell= round(( modangle(ell)/360.)* 1080.)
	ibee= round( (( modangle( bee, 360., negpos=True)+90.)/180.)* 540)

	return tb_408[ibee, iell]

## Convert angle to a specified range by finding the angle modulo the extent of the range. ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def modangle(angle, extent=360., negpos=False):
	offset = 0.
	if(negpos):
		offset = extent/2.

	return ((((angle-offset) % extent) + extent) % extent) - offset


## Get the index of a given velocity ##
 #
 # params list vect           A list
 # params float vel: (One)    Value
 # return int idx             Index of Value in List
 #
 # version 03/2016 
 # author Nguyen Van Hiep ##
def get_index(vect, val):
    idx = (np.abs(np.array(vect)-val)).argmin()
    return idx