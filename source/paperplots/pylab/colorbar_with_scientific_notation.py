import numpy as np
import seaborn as sns
import matplotlib.ticker as tkr
import matplotlib.pyplot as plt
import matplotlib        as mpl

formatter = tkr.ScalarFormatter(useMathText=True)
formatter.set_scientific(True)
formatter.set_powerlimits((-2, 2))

x = np.exp(np.random.uniform(size=(10, 10)) * 10)
sns.heatmap(x, cbar_kws={"format": formatter})
plt.show()


import numpy as np
import pylab as pl

sfmt       = mpl.ticker.ScalarFormatter(useMathText=True) 
# sfmt.set_powerlimits((0, 0))
sfmt.set_scientific(True)
sfmt.set_powerlimits((-2, 2))

x = np.linspace(0,1)
y = np.linspace(0,2)
x,y = np.meshgrid(x,y)
z = 1e-8*(x**2 + y**2)

img1 = pl.contourf(x,y,z)

# cax1    = divider.append_axes("right", size='2%', pad=0.05)

# cbar      =  plt.colorbar(img1, cax=cax1, orientation='vertical')
cbar      =  plt.colorbar(img1, orientation='vertical', format=formatter)
cbar.set_label(label='(10$^{-6}$)',weight='normal')
cbar.ax.tick_params(labelsize=6)

cbar.ax.tick_params(labelsize=12)
cbar.formatter.set_scientific(True)
text = cbar.ax.yaxis.get_offset_text()
text.set_fontsize(14)
text.set_color("blue")
text.set_position((5, 5))

pl.show()