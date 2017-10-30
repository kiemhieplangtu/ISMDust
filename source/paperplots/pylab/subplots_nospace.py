import matplotlib.pyplot as plt

f, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw = {'wspace':0, 'hspace':0})

ax1.grid('on', linestyle='--')
ax1.set_xticklabels([])
ax1.set_yticklabels([])

ax2.grid('on', linestyle='--')
ax2.set_xticklabels([])
ax2.set_yticklabels([])

plt.show()