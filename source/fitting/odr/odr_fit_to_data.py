"""
odr_fit_to_data.py
    A simple program to fit x,y data with uncertainties in both x and y

    Uses Orthogonal Distance Regression

    The example provided is a fit of Gaussian, Lorentzian, or linear functions
    to a data file gauss_xy_error.txt, gauss_dx.txt, or linear_xy_error.txt.

    To fit your own data, you need to change:
    (1) func(p,x) to return the function you are trying to fit,
            p is the parameter vector, x are the independent variable(s)
            Caution: scipy.odr and scipy.optimize.curve_fit require x & p
                in opposite orders.
    (2) the name of the data file read in by numpy.loadtxt,
    (3) the initial guesses (p_guess) for the fit parameters. The fit will not
            converge if your guess is poor.

    For a more accurate estimate of the Cumulative Distribution Function (CDF)
        by Monte Carlo simulation based on the fitted values, set
            Run_Monte_Carlo_CDF_Estimator  = True
            Number_of_MC_iterations        = 1000
        But note that this can take some time (e.g. > 1s per iteration) if the
            x uncertainties are large.

    For information on scipy.odr, see
    http://docs.scipy.org/doc/scipy/reference/odr.html

    Copyright (c) 2011-206 University of Toronto
    Modifications:
        7 February 2016 by David Bailey
            Made notation and uncertainties compatible with curve_fit_to_data,
            in particular prevented uncertainty estimates being reduced if fit
            is "too good'
        7 January 2016 by David Bailey
            Made printing python 3 compatible
        11 January 2013 by David Bailey
    Original Version   :   11 September 2011 by Michael Luzi
    Contact: David Bailey <dbailey@physics.utoronto.ca>
                (www.physics.utoronto.ca/~dbailey)
    License: Released under the MIT License; the full terms are appended
                to the end of this file, and are also available
                at www.opensource.org/licenses/mit-license.php
"""

# Use python 3 consistent printing and unicode
from __future__ import print_function
from __future__ import unicode_literals

import inspect                # http://docs.python.org/library/inspect.html
from matplotlib import pyplot # http://matplotlib.org/api/pyplot_api.html
import numpy                  # http://numpy.scipy.org/
import scipy                  # http://scipy.org/
import scipy.odr, scipy.special, scipy.stats

# Decide if Monte Carlo CDF values are desired
Run_Monte_Carlo_CDF_Estimator = False
Number_of_MC_iterations       = 11

# You can ignore x uncertainties if desired
x_uncertainties_ignored       = False

## define some possible functions to be fitted
##  You may, of course, define your own.
def Cauchy(p,x) :
    # A cauchy / lorentzian / breit wigner peak with contant background:
    B           = p[0]    #   Constant Background
    x1_peak     = p[1]    #   Height above/below background
    x1          = p[2]    #   Central value
    width_x1    = p[3]    #   Half-Width at Half Maximum
    return ( B + x1_peak*(1/(1+((x-x1)/width_x1)**2)) )
def Gauss(p,x) :
    # A gaussian or Normal peak with constant background:
    B           = p[0]    #   Constant Background
    x1_peak     = p[1]    #   Height above/below background
    x1          = p[2]    #   Central value
    width_x1    = p[3]    #   Standard Deviation
    return ( B + x1_peak*numpy.exp(-(x-x1)**2/(2*width_x1**2)) )
def Linear(p,x) :
    # A linear function with:
    #   Constant Background          : p[0]
    #   Slope                        : p[1]
    return p[0]+p[1]*x

## Choose function to fit:
func = Linear
if func.__name__ == "Gauss" :
    # Initial guesses for fit parameters
    p_guess      = (10.,200.,1173.9,0.4,)
    # Gaussian data with neglible x uncertainties
    data_file    = "gauss_dx.txt"
elif func.__name__ == "Cauchy" :
    p_guess      = (10.,200.,1173.9,0.4,)
    # Gaussian data with large x uncertainties
    data_file    = "gauss_xy_errors_ugly.txt"
else :
    # default is linear
    func         = Linear
    p_guess      = (2,2)
    p_guess      = (6.81,0.658)
    # Gaussian data with large x uncertainties
    data_file    = "linear_xy_errors1.txt"
    data_file    = "linear_xy_errors.txt"

## Load data to fit
#     Any lines in the data file starting with '#' are ignored
#     "data" are values, "sigma" are uncertainties
x_data, x_sigma, y_data, y_sigma = numpy.loadtxt(data_file, comments='#', unpack=True, skiprows=0, usecols=None)
if x_uncertainties_ignored:
    # x uncertainties cannot be set exactly to zero with crashing the fit,
    #   but a tiny value seems to do the job.
    x_sigma = len(x_data)*[1e-99]

# Load data for ODR fit
data  = scipy.odr.RealData(x=x_data, y=y_data, sx=x_sigma, sy=y_sigma)
# Load model for ODR fit
model = scipy.odr.Model(func)

## Now fit model to data
#	job=10 selects central finite differences instead of forward
#            differences when doing numerical differentiation of function
#	maxit is maximum number of iterations
#   The scipy.odr documentation is a bit murky, so to be clear note that
#        calling ODR automatically sets full_output=1 in the low level
#        odr function.
fit    = scipy.odr.ODR(data, model, beta0=p_guess, maxit=5000,job=10)

output = fit.run()
p      = output.beta 	   # 'beta' is an array of the parameter estimates
cov    = output.cov_beta   # parameter covariance matrix
#"Quasi-chi-squared" is defined to be the [total weighted sum of squares]/dof
#	i.e. same as numpy.sum((residual/adjusted_err)**2)/dof or
#       numpy.sum(((output.xplus-x)/x_sigma)**2
#                                +((y_data-output.y)/y_sigma)**2)/dof
#	Converges to the conventional chi-squared for zero x uncertainties.
quasi_chisq = output.res_var
uncertainty = output.sd_beta # parameter standard uncertainties
#  scipy.odr appears to scale the parameter uncertainties by quasi_chisq.
#    If the fit is poor, i.e. quasi_chisq is large, the uncertainties
#    are scaled up, which makes sense. If the fit is too good,
#    i.e. quasi_chisq << 1, it suggests that the uncertainties have been
#    overestimated, but it seems risky to scale down the uncertainties. (i.e.
#    The uncertainties shouldn't be too sensistive to random fluctuations of
#    the data, and if the uncertainties are overestimated, this should be
#    understood and corrected
if quasi_chisq < 1.0 :
    uncertainty = uncertainty/numpy.sqrt(quasi_chisq)
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
print("ODR algorithm stop reason: " + output.stopreason[0])
print("\nFit {0} Data points from file: {1}".format(len(x_data),data_file))
print("To Model :")
print(" ",inspect.getsource(func))
if x_uncertainties_ignored:
    print("** WARNING: x uncertainties set to zero in fit. **\n")

print("Estimated parameters and uncertainties")
for i in range(len(p)) :
    print(("   p[{0}] = {1:10.5g} +/- {2:10.5g}"+
           "          (Starting guess: {3:10.5g})").\
            format(i,p[i],uncertainty[i],p_guess[i]))

# number of degrees of freedom for fit
dof = len(x_data) - len(p_guess)

# Note that odr, unlike curve_fit, returns the standard covariance matrix.
print("\nStandard Covariance Matrix : \n", cov, "\n")

print("\nCorrelation Matrix :")
for i,row in enumerate(cov):
    for j in range(len(p)) :
        if (i != j):
            CorrCoeff = cov[i,j]/numpy.sqrt(cov[i,i]*cov[j,j])
        print("{0:< 8.3g}".format(cov[i,j]/numpy.sqrt(cov[i,i]*cov[j,j])),
                end="")
            # Newbie Notes: "{0:< 8.3g}" left justifies output with space in
            #   front of positive numbers, with 3 sig figs;
            #   end="" suppresses new line.
    print()
print("\nCorrelation Coeff. of Parameters (m,b) : \n", CorrCoeff, "\n")
# Calculate initial residuals and the 'adjusted error' for each data point
delta = output.delta  # estimated x-component of the residuals
eps   = output.eps    # estimated y-component of the residuals
# (xstar,ystar) is the point where the 'residual line' (in black)
#   intersects the 'ellipse' created by xerr & yerr.
xstar        = ( x_sigma*numpy.sqrt( ((y_sigma*delta)**2) /
                ( (y_sigma*delta)**2 + (x_sigma*eps)**2 ) ) )
ystar        = ( y_sigma*numpy.sqrt( ((x_sigma*eps)**2) /
                ( (y_sigma*delta)**2 + (x_sigma*eps)**2 ) ) )
adjusted_err = numpy.sqrt(xstar**2 + ystar**2)
# residual is positive if the point lies above the fitted curve,
#             negative if below
residual     = numpy.sign(y_data-func(p,x_data))*numpy.sqrt(delta**2 + eps**2)

print("\nQuasi Chi-Squared/dof   = {0:10.5f}, Chi-Squared CDF = {1:10.5f}%".\
format(quasi_chisq, 100.*float(scipy.special.chdtrc(dof,dof*quasi_chisq))))
print("   WARNING:Above CDF is not valid for large x uncertainties!")

print("\nTo use Monte Carlo simulation to more accurately estimate CDF for")
print('      large x uncertainties, re-run program with ')
print('     "Run_Monte_Carlo_CDF_Estimator = True" and')
print('     "Number_of_MC_iterations >= 1000." This may take some time\n')

print("Run_Monte_Carlo_CDF_Estimator = {0}"\
                  .format(Run_Monte_Carlo_CDF_Estimator))

if Run_Monte_Carlo_CDF_Estimator :
    print("\n**** Running Monte Carlo CDF Estimator ****")
    print("Number_of_MC_iterations = {0}".format(Number_of_MC_iterations),
                end="")
    # Initialize Monte Carlo output distributions
    x_dist           = Number_of_MC_iterations*[None]
    y_dist           = Number_of_MC_iterations*[None]
    p_dist           = Number_of_MC_iterations*[None]
    quasi_chisq_dist = Number_of_MC_iterations*[None]
    # Initialize timing measurement
    import time
    start_time = time.clock()
    for i in range(Number_of_MC_iterations) :
        # Starting with the x and x uncertainty (x_sigma) values from the
        #    data, calculate Monte Carlo values assuming a normal gaussian
        #    distibution.
        x_dist[i]   = numpy.random.normal(x,x_sigma)
        # Calculate y using the Monte Carlo x, then smear by the y uncertainty
        y_dist[i]   = numpy.random.normal(func(p,x),y_sigma)
        # Fit the Monte Carlo x,y pseudo-data to the original model
        data_dist   = scipy.odr.RealData( x=x_dist[i], y=y_dist[i],
                                        sx=x_sigma,  sy=y_sigma)
        model_dist  = scipy.odr.Model(func)
        fit_dist    = scipy.odr.ODR(data_dist, model_dist, p, maxit=5000, job=10)
        output_dist = fit_dist.run()
        p_dist[i]   = output_dist.beta
        quasi_chisq_dist[i] = output_dist.res_var
    end_time = time.clock()

    print(" simulations in {0} seconds.".\
          format(end_time-start_time))
    # Sort the simulated quasi-chi-squared values
    quasi_chisq_sort = numpy.sort(quasi_chisq_dist)
    # Find the index of the sorted simulated quasi-chi-squared value that is
    #	nearest to the data quasi-chi-squared value, to estimate the cdf
    print("\nFraction of quasi-chisquared values larger than observed value:")
    print("    Monte Carlo CDF = {0:6.1f}%".format(
                100.*(1.0-numpy.abs(quasi_chisq_sort-quasi_chisq).argmin()
                                   /float(Number_of_MC_iterations))))
    print("    Minimum, Mean, and Maximum Quasi Chi-Squared values:",end="")
    print("{0:6.2g} {1:6.2g} {2:6.2g}.".\
          format(numpy.min(quasi_chisq_dist),numpy.average(quasi_chisq_dist),
                         numpy.max(quasi_chisq_dist)))
    print("\nAverage and Standard Deviation of MC Fit parameters")
    for i in range(len(p)) :
        print ("   p[{0}] = {1:12.6g} +/- {2:12.6g} ; "+
               "(Min = {3:12.6f}, Max = {4:12.6f}")\
              .format(i,
                      numpy.average(list(zip(*p_dist)[i])),
                      numpy.std(list(zip(*p_dist)[i])),
                      numpy.min(list(zip(*p_dist)[i])),
                      numpy.max(list(zip(*p_dist)[i])) )
    print("\nCheck for any Fit Biases in MC Fit parameters"
            +"(Significance and ratio)")
    # Check if fit on Monte Carlo data tends to give an average output value
    #    for the parameters different from the input parameters.
    for i in range(len(p)) :
        print("   p[{0}] = {1:< 6.2f}   ({2:<12.7g}+/-{3:<12.7g})"\
              .format(i, (numpy.average(list(zip(*p_dist)[i]))/p[i]-1)/
                      (numpy.std(list(zip(*p_dist)[i]))/p[i]/\
                                      numpy.sqrt(Number_of_MC_iterations-1)),
                      numpy.average(list(zip(*p_dist)[i]))/p[i],
                      numpy.std(list(zip(*p_dist)[i]))/p[i]/\
                                      numpy.sqrt(Number_of_MC_iterations-1)))

## Plot

# create figure with light gray background
fig = pyplot.figure(facecolor="0.98")

# 3 rows, 1 column, subplot 1
#   3 rows are declared, but there are only 2 plots; this leaves room for text
#       in the empty 3rd row
fit = fig.add_subplot(211)
# remove tick labels from upper plot (for clean look)
fit.set_xticklabels( () )
pyplot.ylabel("Counts")

pyplot.title("Orthogonal Distance Regression Fit to Data")
# Plot data as red circles, and fitted function as (default) line.
#   For a smooth look,generate many x values for plotting the model
x_model = numpy.arange( min(x_data),max(x_data),
                        (max(x_data)-min(x_data))/1000.)
fit.plot(x_data,y_data,'ro', x_model, func(p,x_model))
# Add error bars on data as red crosses.
fit.errorbar(x_data, y_data, xerr=x_sigma, yerr=y_sigma, fmt='r+')
fit.set_yscale('linear')
#   draw starting guess as dashed green line ('r-')
fit.plot(x_model, func(p_guess,x_model), 'g-', label="Start", linestyle="--")

a = numpy.array([output.xplus,x_data])   # output.xplus = x + delta
b = numpy.array([output.y,y_data])  # output.y = f(p, xfit), or y + epsilon
fit.plot(numpy.array([a[0][0],a[1][0]]), numpy.array([b[0][0],b[1][0]]),
         'k-', label = 'Residuals')
for i in range(1,len(y_data)):
    fit.plot(numpy.array([a[0][i],a[1][i]]), numpy.array([b[0][i],b[1][i]]),
         'k-')
fit.legend(loc='upper left')
fit.grid()

# separate plot to show residuals
residuals = fig.add_subplot(212) # 3 rows, 1 column, subplot 2
residuals.errorbar(x=x_data,y=residual,yerr=adjusted_err,
                   			fmt="r+", label = "Residuals")
# make sure residual plot has same x axis as fit plot
residuals.set_xlim(fit.get_xlim())
# Draw a horizontal line at zero on residuals plot
pyplot.axhline(y=0, color='b')
# Label axes
pyplot.xlabel("energy [KeV]")
# These data look better if 'plain', not scientific, notation is used, and if
#   the tick labels are not offset by a constant (as is done by default).
pyplot.ticklabel_format(style='plain', useOffset=False, axis='x')
pyplot.ylabel("Residuals")

residuals.grid()

pyplot.show()

##### END of odr_fit_to_data.py

"""
Full text of MIT License:

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
"""