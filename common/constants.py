import os
from   basicImport import *

## Angles
DEG2RAD = np.pi/180.

# Some constants
FUKUI   = 4.77e-7 ## tau353 = 4.77e-27 N(HI) 
PLC4PI  = 8.40e-7 ## tau353 = 8.4e-27 N(H) -> Whole sky
PLCLOW  = 6.60e-7 ## tau353 = 6.6e-27 N(H) -> Low NHI

## Paths
HOME    = os.getenv("HOME")
APPPATH = HOME    + '/PhD@MQ/projects/'
HJ1PATH = APPPATH + 'ISMDust/'
HJ2PATH = APPPATH + 'dark/'
DATPATH = HOME    + '/hdata/'

## Scaling HT03 factor
HT03SCL = 1.26 #1.26 1.14
NHICST  = 1.93988 # NHI/e18 = 1.93988 tau*Ts*Wid