## Basic packages
from   basicImport             import *

from   astropy.io              import fits
from   mpl_toolkits.axes_grid1 import make_axes_locatable ## for Plotting L.O.S Samples 

## Add modules
from   restore                 import restore
import imgData                 as     imgdat
import linFit                  as     lnfit
import txtData                 as     txtDat
import statistics              as     stats
import hpxMap                  as     hpx

## Add style ##
 # 1. Check dir: matplotlib.get_configdir() ; it's usually: '/home/vnguyen/.config/matplotlib'
 # 2. mkdir stylelib
 # 3. create files: stylelib/mystyle.mplstyle
plt.style.use('mystyle')