# Use python 3 consistent printing and unicode
from __future__ import print_function
from __future__ import unicode_literals

from basicImport import *
from restore     import restore
from scipy.odr   import *


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
		parinfo.append(cp.deepcopy(parbase))

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
		parinfo.append(cp.deepcopy(parbase))

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