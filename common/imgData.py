from basicImport import *

#######################################################################
# End - Block for List
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################

## Find common elements of 2 lists ##
 #
 # params list x 
 # params list y 
 #
 # return list List of common elements
 #
 # version 12/2017 
 # author Nguyen Van Hiep ##
def common_elements(x, y):
 	return list(set(x) & set(y))

## Find Different elements of 2 lists ##
 #
 # params list x 
 # params list y 
 #
 # return list List of diff elements
 #
 # version 12/2017 
 # author Nguyen Van Hiep ##
def diff_elements(x, y):
 	return list(set(x) - set(y))

#######################################################################
# End - Block for List
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################







#######################################################################
# Block for Writing to file
# Author: Van Hiep
# Version: 08/2017
##################### - START - #######################################

## Write to file ##
 #
 # params str filename
 # params str string to write to file
 #
 # Example of a string ##
 # string = str(n)+ '\t'+ sc+ '\t'+ str( round(noh[sc]['l'],8) )+ '\t' + '\n'
 #
 # return Void
 #
 # version 08/2017 
 # author Nguyen Van Hiep ##
def write2file(filename, string):
	file = open(filename,'w') 
	file.write( string )
	file.close()

## Write to screen ##
 #
 # params str  style  Format to write
 # params list cols   Columns to write
 #
 # Example of a string ##
 # print '{:10s}{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}\t{:08.4f}'\
 # 	.format(src, err, nh_er, nh2_er, nhier, n_h, nhi, nh2, nh_er, nh_er, nhier, nh2_er)
 #
 # return Void
 #
 # version 08/2017 
 # author Nguyen Van Hiep ##
def printCols(style='{:10s}{:08.4f}\t{:08.4f}', cols=[]):
	print (style)
	print (cols)
	# print '\''+style+'\''.format(cols)

#######################################################################
# End - Block for Writing to file
# Author: Van Hiep
# Version: 08/2017
##################### - END - #########################################


## Get x,y ID ranges for WCO contours #
 #
 # params float x1 X-min
 # params float x2 X-max
 # params float y1 Y-min
 # params float y2 Y-max
 #
 # return list of ID x-y ranges
 # 
 # version 9/2017
 # Author Van Hiep ##
def get_xy_id_ranges(x1=-180.,x2=180.,y1=-80.,y2=89.):
	x = np.arange(180.,-180.125,-0.125)
	y = np.arange(-80.,89.125,0.125)

	# xid1, = np.where(x==x2)
	# xid2, = np.where(x==x1)
	# yid1, = np.where(y==y1)
	# yid2, = np.where(y==y2)

	xid1 = (np.abs(x-x2)).argmin()
	xid2 = (np.abs(x-x1)).argmin()
	yid1 = (np.abs(y-y1)).argmin()
	yid2 = (np.abs(y-y2)).argmin()

	print('x1,x2,y1,y2', x1,x2,y1,y2)

	print('xid1, xid2, yid1, yid2', xid1, xid2, yid1, yid2)

	return yid1, yid2, xid1, xid2

## WCO map from Dame #
 #
 # params 
 # return 2-D array X, Y, Z
 # 
 # version 9/2017
 # Author Van Hiep ##
def Wco_map():
	from   astropy.io import fits

	image_file = os.getenv("HOME") + '/hdata/co/mlat_Wco_mom.fits'
	hdulist    = fits.open(image_file)

	# print hdulist[0].header
	data       = hdulist[0].data[:, ::-1]
	# data     = hdulist[0].data
	print('WCO data shape', data.shape )

	aa   = data[:,960:2881]
	bb   = data[:, 0:960]
	cc   = np.concatenate((aa, bb), axis=1)
	data = cc	


	Nx   = 2881
	Ny   = 1353

	x    = np.arange(180.,-180.125,-0.125)
	y    = np.arange(-80.,89.125,0.125)
	X, Y = np.meshgrid(x, y)

	return X,Y,data

## Patch of WCO map #
 #
 # params 2D-Arr X X-array
 # params 2D-Arr Y X-array
 # params 2D-Arr Z Z-array
 # params float x1 x-min
 # params float x2 x-max
 # params float y1 y-min
 # params float y2 y-max
 #
 # return list of ID x-y ranges
 # 
 # version 9/2017
 # Author Van Hiep ##
def patch_Wco_map(X,Y,Z, xmin=-180., xmax=180., ymin=-80., ymax=89.):
	# X,Y,Z       = Wco_map()
	x1,x2,y1,y2 = get_xy_id_ranges(xmin, xmax, ymin, ymax)

	xd = X[x1:x2,y1:y2]
	yd = Y[x1:x2,y1:y2]
	zd = Z[x1:x2,y1:y2]

	return xd, yd, zd


## Retreive a SINGLE value of 408 t_b from haslam et al. ##
 #
 # params float ell Galactic-longitude
 # params float bee Galactic-latitude
 #
 # return float Tbg_408 Background-Temperature at (l,b)
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def get_tb_408(ell,bee,tb_408):
	iell= round(( modangle(ell)/360.)* 1080.)
	ibee= round( (( modangle( bee, 360., negpos=True)+90.)/180.)* 540)

	return tb_408[ibee, iell]

## Convert angle to a specified range by finding the angle modulo the extent of the range. ##
 #
 # params 
 # params 
 #
 # return 
 #
 # version 08/2016 
 # author Nguyen Van Hiep ##
def modangle(angle, extent=360., negpos=False):
	offset = 0.
	if(negpos):
		offset = extent/2.

	return ((((angle-offset) % extent) + extent) % extent) - offset


## Get the index of a given velocity ##
 #
 # params list vect           A list
 # params float vel: (One)    Value
 # return int idx             Index of Value in List
 #
 # version 03/2016 
 # author Nguyen Van Hiep ##
def get_index(vect, val):
    idx = (np.abs(np.array(vect)-val)).argmin()
    return idx