import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md
import copy

from numpy               import array
from restore             import restore
from plotting            import cplot
from mpfit               import mpfit
from scipy.odr           import *


x = [3., 0., 3., 1.5, -4., -3., 3., 4.2, 2.7]
y = [4., -1., 4., 1., -4., -4., 2.8, 4., 2.5]

coxy = md.cov_xy(x,y)
print coxy
varx = md.var(x)
vary = md.var(y)

rho = coxy/varx/vary
print rho

plt.plot(np.array(x)/np.array(varx), np.array(y)/np.array(vary), 'r*')
plt.xlabel('X/varx')
plt.ylabel('Y/vary')
plt.show()