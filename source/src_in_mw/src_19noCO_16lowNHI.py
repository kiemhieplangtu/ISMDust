import sys, os
sys.path.insert(0, os.getenv("HOME")+'/dark/common') # add folder of Class
import numpy             as np
import healpy            as hp
import matplotlib.pyplot as plt
import module            as md

def plot_mwd(org=0,title='', projection='mollweide'):
    ''' RA, Dec are arrays of the same length.
    RA takes values in [0,360), Dec in [-90,90],
    which represent angles in degrees.
    org is the origin of the plot, 0 or a multiple of 30 degrees in [0,360).
    title is the title of the figure.
    projection is the kind of projection: 'mollweide', 'aitoff', 'hammer', 'lambert'
    '''
    info         = md.read_19src_noco_nooh(fname = '../oh/result/19src_noCO_noOH.txt')
    lownhi       = md.read_23rc_lownhi(fname = '../oh/result/16src_lowNHI.txt')
    xtick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    ytick_labels = np.array([-75, -60, -45, -30, -15, 0, 15, 30, 45, 60, 75])
    # tick_labels = np.remainder(tick_labels+360+org,360)


    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(111, projection=projection, axisbg ='LightCyan')
    # ax.scatter(np.radians(x),np.radians(Dec))  # convert degrees to radians
    for i in range(0, len(lownhi['src'])):
		src   = lownhi['src'][i]
		l     = lownhi['l'][i]
		if (l<180.):
			l = -l
		else:
			l = -l + 360.

		b     = lownhi['b'][i]
		print 'aa', l
		ax.scatter(np.radians([l]),np.radians([b]), color='r', marker='^', s=90)
		ax.text(np.radians([l]),np.radians([b]), '  '+src, fontsize=8, position=(np.radians([l+1.]),np.radians([b-2.])) )

    for i in range(0, len(info['src'])):
		src   = info['src'][i]
		l     = info['l'][i]
		if (l<180.):
			l = -l
		else:
			l = -l + 360.

		b     = info['b'][i]
		print 'bb', l
		ax.scatter(np.radians([l]),np.radians([b]), marker='h', s=90)
		ax.text(np.radians([l]),np.radians([b]), src, fontsize=8)

    ax.set_xticklabels(xtick_labels, fontweight='bold')     # we add the scale on the x axis
    ax.set_yticklabels(ytick_labels, fontweight='bold') 
    ax.set_title(title)
    ax.title.set_fontsize(15)
    ax.set_xlabel("Galactic longitude ($^{o}$)")
    ax.xaxis.label.set_fontsize(24)
    ax.xaxis.label.set_fontweight('bold')
    ax.set_ylabel("Galactic latitude ($^{o}$)")
    ax.yaxis.label.set_fontsize(24)
    ax.yaxis.label.set_fontweight('bold')
    ax.grid(True)

## ============== ##
info   = md.read_19src_noco_nooh(fname = '../oh/result/19src_noCO_noOH.txt')
lownhi = md.read_23rc_lownhi(fname = '../oh/result/16src_lowNHI.txt')

plot_mwd()

plt.tight_layout()
plt.savefig('src_locations.eps', format='eps', dpi=600)
plt.show()