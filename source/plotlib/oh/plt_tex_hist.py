import sys, os
sys.path.insert(0, os.getenv("HOME")+'/ISMDust/common') # add folder of Class

import numpy             as np
import matplotlib.pyplot as plt
import matplotlib
import pandas            as pd

from scipy.optimize import curve_fit
from scipy.io.idl   import readsav
from numpy          import array
from restore        import restore
from plotting       import cplot
from gauss_fit      import gfit

## Create a line ##
 #
 # params list x x-data
 # params list y y-data
 # params string label Label of line
 # params dict prop Properties of line
 # return dict ret All infor of line
 #
 # version 07/2016 
 # author Nguyen Van Hiep ##
def daddsada(tb):
	return ''

#============== MAIN ==============#
# dat                 = readsav('data/doall_fitgauss_stokesiover2.sav')
dat                 = readsav('data/doall_fitgauss.sav')

## For OH1665 Absorption 
src1                = dat.amg_1665.sname
ell1                = dat.amg_1665.ell
bee1                = dat.amg_1665.bee
ng1                 = dat.amg_1665.nrgauss
zro1                = dat.amg_1665.zrocnm
taucnm1             = dat.amg_1665.taucnm
cencnm1             = dat.amg_1665.cencnm
widcnm1             = dat.amg_1665.widcnm
tspincnm1           = dat.amg_1665.tspincnm
sigzrocnm1          = dat.amg_1665.sigzrocnm
sigtaucnm1          = dat.amg_1665.sigtaucnm
sigcencnm1          = dat.amg_1665.sigcencnm
sigwidcnm1          = dat.amg_1665.sigwidcnm
sigtspincnm1        = dat.amg_1665.sigtspincnm
continuum_defln1    = dat.amg_1665.continuum_defln
sigcontinuum_defln1 = dat.amg_1665.sigcontinuum_defln
continuum_em1       = dat.amg_1665.continuum_em


## For OH1667 Absorption 
src2                = dat.amg_1667.sname
ng2                 = dat.amg_1667.nrgauss
ell2                = dat.amg_1667.ell
bee2                = dat.amg_1667.bee
zro2                = dat.amg_1667.zrocnm
taucnm2             = dat.amg_1667.taucnm
cencnm2             = dat.amg_1667.cencnm
widcnm2             = dat.amg_1667.widcnm
tspincnm2           = dat.amg_1667.tspincnm
sigzrocnm2          = dat.amg_1667.sigzrocnm
sigtaucnm2          = dat.amg_1667.sigtaucnm
sigcencnm2          = dat.amg_1667.sigcencnm
sigwidcnm2          = dat.amg_1667.sigwidcnm
sigtspincnm2        = dat.amg_1667.sigtspincnm
continuum_defln2    = dat.amg_1667.continuum_defln
sigcontinuum_defln2 = dat.amg_1667.sigcontinuum_defln
continuum_em2       = dat.amg_1667.continuum_em

## For OH1665 Emmission
esrc1                = dat.emg_1665.sname
# ell2                 = dat.emg_1665.ell
# bee2                 = dat.emg_1665.bee
eng1                 = dat.emg_1665.nrgauss
bsl1                 = dat.emg_1665.tbaseline
cont_em1             = dat.emg_1665.continuum_em
hgtwnm1              = dat.emg_1665.hgtwnm
cenwnm1              = dat.emg_1665.cenwnm
widwnm1              = dat.emg_1665.widwnm
fwnm1                = dat.emg_1665.fwnm
sigcont1             = dat.emg_1665.sigcontinuum
sighgtwnm1           = dat.emg_1665.sighgtwnm
sigcenwnm1           = dat.emg_1665.sigcenwnm
sigwidwnm1           = dat.emg_1665.sigwidwnm
sigfwnm1             = dat.emg_1665.sigfwnm

## For OH1667 Emmission 
esrc2                = dat.emg_1667.sname
# ell2                 = dat.emg_1667.ell
# bee2                 = dat.emg_1667.bee
eng2                 = dat.emg_1667.nrgauss
bsl2                 = dat.emg_1667.tbaseline
cont_em2             = dat.emg_1667.continuum_em
hgtwnm2              = dat.emg_1667.hgtwnm
cenwnm2              = dat.emg_1667.cenwnm
widwnm2              = dat.emg_1667.widwnm
fwnm2                = dat.emg_1667.fwnm
sigcont2             = dat.emg_1667.sigcontinuum
sighgtwnm2           = dat.emg_1667.sighgtwnm
sigcenwnm2           = dat.emg_1667.sigcenwnm
sigwidwnm2           = dat.emg_1667.sigwidwnm
sigfwnm2             = dat.emg_1667.sigfwnm

for i in range(len(src1)):
	# print src1[i], taucnm1[i], cencnm1[i], widcnm1[i], tspincnm1[i]
	print i, src1[i], continuum_em1[i], continuum_em2[i], continuum_defln1[i], continuum_defln2[i], cont_em1[i], cont_em2[i]

sys.exit()

## Plot

# matplotlib.style.use('grayscale')
tbg = list(set(continuum_em1))
tbg = [3.41378, 4.51589,3.23905,3.66018,3.44514,3.70722, 3.6669,3.54818, 3.53922,3.57058]

print 'mean Tex1, mean Tex2'
print np.mean( np.array(tspincnm1) )
print np.mean( np.array(tspincnm2) )
print len(tbg)
print 'Rewrite Tbg:'
print tbg

## Overlapping
plt.rc('font', weight='bold')
plt.rc('text', usetex=True)
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sfmath} \boldmath']

c1           = 'r'
c2           = 'dimgrey'
c3           = 'black'
fig          = plt.figure(figsize=(8,11.2))
ax           = fig.add_subplot(111); ax.set_rasterized(True)
major_ticks  = np.arange(0, 100, 5)                                              
minor_ticks  = np.arange(0, 100, 1)
major_yticks = np.arange(0., 25, 2)                                              
minor_yticks = np.arange(0., 25, 1)
plt.hist(tspincnm1, alpha=0.4,  label=r'${\rm{T_{ex}(1665)}}$', color=c1, ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,26.0), bins=13, lw=3, edgecolor=c1)
plt.hist(tspincnm2, alpha=0.3, label=r'${\rm{T_{ex}(1667)}}$', color=c2, ls='-',  histtype='stepfilled', stacked=False, fill=True, range=(0.0,26.0), bins=13, lw=3, edgecolor=c3)
plt.hist(tbg,       alpha=1.0,  label=r'${\rm{T_{bg}}}$',       color=c3, ls='--', histtype='step',       stacked=False, fill=False, range=(0.0,26.0), bins=13, lw=3, edgecolor=c3)

plt.hist(tspincnm1, color=c1, ls='-',  histtype='stepfilled', stacked=False, fill=False, range=(0.0,26.0), bins=13, lw=3, edgecolor=c1)
plt.hist(tspincnm2, color=c2, ls='-',  histtype='stepfilled', stacked=False, fill=False, range=(0.0,26.0), bins=13, lw=3, edgecolor=c3)
plt.hist(tbg,       color=c3, ls='--', histtype='step',       stacked=False, fill=False, range=(0.0,26.0), bins=13, lw=4, edgecolor=c3)

bins             = np.arange(0.,26.1,2.)
tex1, bin_edges1 = np.histogram(tspincnm1, bins=bins, range=(0.0,26.0), normed=False)
tex2, bin_edges2 = np.histogram(tspincnm2, bins=bins, range=(0.0,26.0), normed=False)
bins             = np.arange(0.,26.1,2.)
tbg,  bin_edges3 = np.histogram(tbg,       bins=bins, range=(0.0,26.0), normed=False)

bincen1 = 0.5*(bin_edges1[1:] + bin_edges1[:-1])
bincen2 = 0.5*(bin_edges2[1:] + bin_edges2[:-1])
bincen3 = 0.5*(bin_edges3[1:] + bin_edges3[:-1])

# plot1   = ax.plot(bincen1, tex1, ls='',  color=c1, marker='d',  markersize=12)
# plot2   = ax.plot(bincen2, tex2, ls='',  color=c2, marker='d',  markersize=12)
# plot3   = ax.plot(bincen3, tbg,  ls='',  color=c3, marker='d',  markersize=12)

plot1   = ax.plot(bincen1, tex1, ls='',  color=c1, marker='*',  markersize=16)
plot2   = ax.plot(bincen2, tex2, ls='',  color=c2, marker='d',  markersize=14)
plot3   = ax.plot(bincen3, tbg,  ls='',  color=c3, marker='^',  markersize=14)

plt.title('', fontsize=22)
plt.ylabel(r'\textbf{Number}', fontsize=22,fontweight='bold')
plt.xlabel(r'\textbf{Temperature (K)}', fontsize=22, fontweight='bold')
ax.set_xticks(major_ticks)                                                       
ax.set_xticks(minor_ticks, minor=True)                                           
ax.set_yticks(major_yticks)                                                       
ax.set_yticks(minor_yticks, minor=True)
plt.tick_params(axis='x', labelsize=22, pad=8)
plt.tick_params(axis='y', labelsize=22)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=22)
plt.grid(False)
plt.xlim(0.,26.)
plt.ylim(0.,12.)

# plt.tight_layout()

plt.savefig('Tex_hist.eps', format='eps', dpi=600)
plt.show()




sys.exit()
## Use numpy histogram, no-color
tbg              = list(set(continuum_em1))
fig              = plt.figure(figsize=(12,16))
ax               = fig.add_subplot(111)
major_ticks      = np.arange(0, 100, 5)                                              
minor_ticks      = np.arange(0, 100, 1) 

bins             = np.arange(0.,20.1,2.)
tex1, bin_edges1 = np.histogram(tspincnm1, bins=bins, range=(0.0,20.0), normed=False)
tex2, bin_edges2 = np.histogram(tspincnm2, bins=bins, range=(0.0,20.0), normed=False)
bins             = np.arange(0.,26.1,2.)
tbg,  bin_edges3 = np.histogram(tbg,       bins=bins, range=(0.0,26.0), normed=False)

bincen1 = 0.5*(bin_edges1[1:] + bin_edges1[:-1])
bincen2 = 0.5*(bin_edges2[1:] + bin_edges2[:-1])
bincen3 = 0.5*(bin_edges3[1:] + bin_edges3[:-1])

plot1   = ax.plot(bincen1, tex1, ls='-', color='k', label='$T^{OH^{1665}}_{ex}$', marker='s', markersize=8,  lw=2)
plot2   = ax.plot(bincen2, tex2, ls='--',  color='k', label='$T^{OH^{1667}}_{ex}$', marker='*', markersize=13, lw=2)
plot3   = ax.plot(bincen3, tbg,  ls='-.',  color='k', label='$T_{bg}$',             marker='D', markersize=9,  lw=2)

plt.title('', fontsize=30)
plt.ylabel('Numbers', fontsize=35)
plt.xlabel('Temperature (K)', fontsize=35)
ax.set_xticks(major_ticks)                                                       
ax.set_xticks(minor_ticks, minor=True)                                           
ax.set_yticks(major_ticks)                                                       
ax.set_yticks(minor_ticks, minor=True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(0.,26.)
plt.ylim(0.,20.)
plt.savefig('t2.eps', format='eps', dpi=1000)
plt.show()

## Use numpy histogram, with colors
fig              = plt.figure(figsize=(12,16))
ax               = fig.add_subplot(111)
major_ticks      = np.arange(0, 100, 5)                                              
minor_ticks      = np.arange(0, 100, 1)

ax.plot(bincen1, tex1, ls='-', color='k', label='$T^{OH_{ex}(1665)$', marker='s', markersize=8, lw=2)
ax.plot(bincen2, tex2, ls='--',  color='r', label='$T_{ex}(1667)$',   marker='*', markersize=13, lw=2)
ax.plot(bincen3, tbg,  ls='-.',  color='g', label='$T_{bg}$',         marker='D', markersize=9, lw=2)

plt.title('', fontsize=30)
plt.ylabel('Numbers', fontsize=35)
plt.xlabel('Temperature (K)', fontsize=35)
ax.set_xticks(major_ticks)                                                       
ax.set_xticks(minor_ticks, minor=True)                                           
ax.set_yticks(major_ticks)                                                       
ax.set_yticks(minor_ticks, minor=True)
plt.tick_params(axis='x', labelsize=20)
plt.tick_params(axis='y', labelsize=20)
plt.tick_params(which='both', width=2)
plt.tick_params(which='major', length=9)
plt.tick_params(which='minor', length=4)
plt.legend(loc='upper right', fontsize=18)
plt.grid(False)
plt.xlim(0.,26.)
plt.ylim(0.,20.)
plt.savefig('t3.eps', format='eps', dpi=1000)
plt.show()