import math
import sys
import numpy as np

# Arguments of do_block_fes.py
# - FILE: input file, 1 column per CV + weights (optional)
# - NCV: number of CVs
# - *MIN: minimum value of CV
# - *MAX: max value of CV
# - *NBIN: # points in output free energy
# - KBT: temperature in energy units (kJoule/mol)
# - N: Block size
#
# * = repeat this block for each CV
# Example with 2 CVs:
# python3 do_block_fes.py phi_psi_w.dat 2 -3.141593 3.141593 50 -3.141593 3.141593 50 2.494339 100
#
#
# Author: Max Bonomi (mbonomi@pasteur.fr)
#
# 
# useful functions
# nD indexes from 1D index
def get_indexes_from_index(index, nbin):
    indexes = []
    # get first index
    indexes.append(index%nbin[0])
    # loop
    kk = index
    for i in range(1, len(nbin)-1):
        kk = ( kk - indexes[i-1] ) / nbin[i-1]
        indexes.append(kk%nbin[i])
    if(len(nbin)>=2):
      indexes.append( ( kk - indexes[len(nbin)-2] ) / nbin[len(nbin) -2] )
    return tuple(indexes)
 
# nD indexes from values
def get_indexes_from_cvs(cvs, gmin, dx, nbin):
    idx = []
    for i in range(0, len(cvs)):
        j = int( round( ( cvs[i] - gmin[i] ) / dx[i] ) )
        # check boundaries
        if(j>=nbin[i]):
          print("Point outside grid, check boundaries!")
          exit()
        idx.append(j)
    return tuple(idx)

# 1D index from values
def get_index_from_cvs(cvs, gmin, dx, nbin):
    # get nD indices from value 
    idx = get_indexes_from_cvs(cvs, gmin, dx, nbin)
    # transform in 1D index
    i = idx[-1]
    for j in range(len(nbin)-1,0,-1):
        i = i*nbin[j-1]+idx[j-1]
    return i

# grid points from nD indexes
def get_points_from_indexes(idx, gmin, dx):
    xs = []
    for i in range(0, len(idx)):
        xs.append(gmin[i] + float(idx[i]) * dx[i])
    return xs

# read CV/weight file and create numpy arrays
def read_file(filename,gmin,dx,nbin):
    # read file and store lists 
    cvs=[]; ws=[]
    # number of cv
    ncv = len(gmin)
    for lines in open(filename, "r").readlines():
        riga = lines.strip().split()
        # check format
        if(len(riga)!=ncv and len(riga)!=ncv+1):
          print (filename,"is in the wrong format!")
          exit()
        # read CVs
        cv = []
        for i in range(0, ncv): cv.append(float(riga[i]))
        # get index in flattened array
        idx = get_index_from_cvs(cv, gmin, dx, nbin)
        # read weight, if present
        if(len(riga)==ncv+1):
          w = float(riga[ncv])
        else: w = 1.0
        # store into cv and weight lists
        cvs.append(idx)
        ws.append(w)
    # return numpy arrays
    return np.array(cvs),np.array(ws)

# 1) READ INPUT parameters
# FILE with CV trajectory and (optionally) weights
FILENAME_ = sys.argv[1]
# number of CVs 
NCV_ = int(sys.argv[2])
# read minimum, maximum and number of bins for FES grid
gmin = []; gmax = []; nbin = []
for i in range(0, NCV_):
    i0 = 3*i + 3 
    gmin.append(float(sys.argv[i0]))
    gmax.append(float(sys.argv[i0+1]))
    nbin.append(int(sys.argv[i0+2]))
# read KBT_
KBT_ = float(sys.argv[3*NCV_+3])
# block size 
BSIZE_ = int(sys.argv[-1])

# 2) SETUP
# define bin sizes
dx = []
for i in range(0, NCV_):
    dx.append( (gmax[i]-gmin[i])/float(nbin[i]-1) )
# total numbers of bins
nbins = 1
for i in range(0, len(nbin)): nbins *= nbin[i]
# read file and store arrays
cv, w = read_file(FILENAME_, gmin, dx, nbin)
# total number of data points
ndata = cv.shape[0]
# number of blocks
nblock = int(ndata/BSIZE_)
# prepare numpy arrays for histogram and normalization
histo = np.zeros((nbins,nblock))
norm  = np.zeros(nblock)

# 3) FILL IN HISTOGRAM ARRAY
for iblock in range(0, nblock):
    # define range
    i0 = iblock * BSIZE_ 
    i1 = i0 + BSIZE_
    # cycle on points in the block
    for i in range(i0, i1):
        # update histogram
        histo[cv[i],iblock] += w[i]
    # calculate normalization of the block
    norm[iblock] = np.sum(w[i0:i1])
    # normalize block
    histo[:,iblock] /= norm[iblock]

# 4) CALCULATE WEIGHTED AVERAGE AND VARIANCE
# now we calculate weighted average across blocks
ave   = np.sum(histo*norm, axis=1) / np.sum(norm)
avet  = np.transpose(np.tile(ave, (nblock,1)))
# and variance
var = np.sum(np.power( norm * (histo-avet), 2), axis=1) / np.power(np.sum(norm), 2)

# 5) PRINT FES + ERROR
log = open("fes."+str(BSIZE_)+".dat", "w")
# this is needed to add a blank line
xs_old = []
for i in range(0, nbins):
    # get the indexes in the multi-dimensional grid
    idx = get_indexes_from_index(i, nbin)
    # get values for grid point
    xs = get_points_from_indexes(idx, gmin, dx)
    # add a blank line for gnuplot
    if(i == 0):
      xs_old = xs[:] 
    else:
      flag = 0
      for j in range(1,len(xs)):
          if(xs[j] != xs_old[j]):
            flag = 1
            xs_old = xs[:] 
      if (flag == 1): log.write("\n")
    # print grid point 
    for x in xs:
        log.write("%12.6lf " % x)
    # calculate fes and error
    if (ave[i]>0):
       # fes
       fes = -KBT_ * math.log(ave[i])
       # variance fes
       varf = math.pow( KBT_ / ave[i], 2.0) * var[i]
       # error fes 
       errf = math.sqrt(varf)
       # printout
       log.write("   %12.6lf %12.6lf\n" % (fes, errf))
    else:
       log.write("   %12s %12s\n" % ("Inf","Inf"))
log.close()
