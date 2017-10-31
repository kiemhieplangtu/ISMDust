import os, sys, shutil
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import module            as md


from   astropy.io        import fits
from astropy.table       import Table


############## MAIN ###########
X,Y,Z       = md.Wco_map()
xd, yd, zd  = md.patch_Wco_map(X,Y,Z)



# plt.figure(2, figsize=(12,12))
# # CS = plt.contour(X, Y, np.log10(data), 5)
# # CS = plt.contour(np.log10(data), 5)
# CS = plt.contour(xd, yd, np.log10(zd), 5)
# plt.clabel(CS, inline=1, fontsize=10)
# plt.title('WCO12')
# plt.xlim(180., -180.)
# plt.ylim(-80., 89)


plt.figure(2, figsize=(12,12))
# CS = plt.contour(X, Y, np.log10(data), 5)
# CS = plt.contour(np.log10(data), 5)
CS = plt.contour(xd, yd, zd)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('WCO12')

plt.show()