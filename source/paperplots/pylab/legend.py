import matplotlib.pyplot as plt
import numpy as np

subs=['one','two','three']
x=[1,2,3]
y=[1,2,3]
yerr=[2,3,1]
xerr=[0.5,1,1]
fig,(ax1)=plt.subplots(1,1)

for i in np.arange(len(x)):
    ax1.errorbar(x[i],y[i],yerr=yerr[i],xerr=xerr[i],label=subs[i],ecolor='black',marker='o',ls='')

# get handles
handles, labels = ax1.get_legend_handles_labels()
# remove the errorbars
handles = [h[0] for h in handles]
# use them in the legend
ax1.legend(handles, labels, loc='upper left',numpoints=1)


plt.show()