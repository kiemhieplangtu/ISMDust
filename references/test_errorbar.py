import matplotlib.pyplot as plt
import numpy as np

s=[19.0, 20.0, 21.0, 22.0, 24.0]
v=np.array([36.5, 66.814250000000001, 130.17750000000001, 498.57466666666664, 19.41])
verr=np.array([9., 9., 9., 9., 9.])
plt.errorbar(v,s,xerr=verr)
# plt.ylim(1E1,1E4)
plt.yscale('log')
plt.xscale('log')
plt.show()