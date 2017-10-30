import matplotlib.pyplot as plt
import numpy             as np

from numpy  import array
import operator

# Read glon & glat #
#
# params string fname Filename
#
# return dict info of src
# 
# Author Van Hiep
##
def read_lb(fname = '../dust/ebv2nh/data/ebv_sfd98_sf2011.txt'):
	ret = {}

	ret['l0']   = [] # 0 - no CO
	ret['b0']   = []
	ret['l1']   = [] # 1 - Low N(HI)
	ret['b1']   = []

	ret['src']  = []

	file = open (fname,'r')
	file.readline() # comment line
	file.readline() # comment line

	for line in file:
	    line    = line.strip()
	    columns = line.split()

	    l = float(columns[2])
	    if (l > 180.) :
	    	l = l - 360.

	    typ = int(columns[0])
	    if (typ == 0) :
	    	ret['l0'].append(l)
	    	ret['b0'].append(float(columns[3]))
	    elif (typ == 1) :
	    	ret['l1'].append(l)
	    	ret['b1'].append(float(columns[3]))

	    ret['src'].append(columns[1])

	file.close()

	return ret

# plot sorces in MW map #
#
# params dict data data to plot
#
# return void
# 
# Author Van Hiep
##
def map_src_in_mw(data):

	l0 = data['l0'] # 
	b0 = data['b0'] # 26 src without CO

	l1 = data['l1'] # 
	b1 = data['b1'] # Low N(HI)

	t = np.arange(-180., 180.0, 0.1)

	fig = plt.figure(figsize=(18,10))
	ax  = fig.add_subplot(111); ax.set_rasterized(True)

	plt.plot(l0,b0, 'r^', label='without CO', markersize=10)
	plt.plot(l1,b1, 'b^', label='low N(HI)', markersize=10)
	plt.plot(t,0.*t, 'k--', label='', linewidth=0.6)
	plt.plot(t,0.*t+5., 'k-.', label='', linewidth=0.4)
	plt.plot(t,0.*t-5., 'k-.', label='', linewidth=0.4)
	plt.plot(t,0.*t+10., 'k:', label='', linewidth=0.4)
	plt.plot(t,0.*t-10., 'k:', label='', linewidth=0.4)

	plt.xlabel('Galactic longitude $(^{o})$', fontsize=35)
	plt.ylabel('Galactic latitude $(^{o})$', fontsize=35)
	plt.title('Locations of 26 LOS without CO and 23 low N(HI) LOS', fontsize=30)
	plt.xlim(-180., 180.)
	plt.ylim(-90., 90.)
	# plt.grid(True)
	plt.tick_params(axis='x', labelsize=18)
	plt.tick_params(axis='y', labelsize=18)

	# plt.text(0.21, 1.31, 'Arecibo beam at 1.4GHz = 3.5\'', color='blue', fontsize=12)
	# plt.text(0.21, 0.92, 'a = '+str(a)+'  b = '+str(b), color='blue', fontsize=12)

	plt.tight_layout()
	plt.legend(loc='upper right', fontsize=18)
	plt.show()

#================= MAIN ========================#

data = read_lb()
map_src_in_mw(data)