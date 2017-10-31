import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import matplotlib.pyplot as plt
import numpy             as np
import matplotlib        as mpl
import healpy            as hp
import pylab             as pl
import module            as md

from numpy               import array
from restore             import restore
from plotting            import cplot
from mpfit               import mpfit

specs = md.read_hi_specs(fname = 'data/nhi_opac_specs.txt')
v     = specs['3C245']['v']
Texp  = specs['3C245']['Texp']
tau   = specs['3C245']['tau']

plt.plot(v, tau)
plt.show()

plt.plot(v, Texp)
plt.show()