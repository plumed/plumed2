import math
import sys

# read FILE with CVs and weights
FILENAME_ = sys.argv[1]
# number of CVs for FES
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
# read output fes
FESFILE_ = sys.argv[3*NCV_+4]

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
    return indexes 

def get_indexes_from_cvs(cvs, gmin, dx):
    keys = []
    for i in range(0, len(cvs)):
        keys.append(int( round( ( cvs[i] - gmin[i] ) / dx[i] ) ))
    return tuple(keys)

def get_points(key, gmin, dx):
    xs = []
    for i in range(0, len(key)):
        xs.append(gmin[i] + float(key[i]) * dx[i])
    return xs

# define bin size
dx = []
for i in range(0, NCV_):
    dx.append( (gmax[i]-gmin[i])/float(nbin[i]-1) )

# create histogram
histo = {}

# read file and fill in histogram
for lines in open(FILENAME_, "r").readlines():
    riga = lines.strip().split()
    # check format
    if(len(riga)!=NCV_ and len(riga)!=NCV_+1):
      print FILENAME_,"is in the wrong format!"
      exit()
    # read CVs
    cvs = []
    for i in range(0, NCV_):
        cvs.append(float(riga[i]))
    # get indexes
    key = get_indexes_from_cvs(cvs, gmin, dx)
    if(len(riga)==NCV_+1):
      # read weight
      w = float(riga[NCV_])
    else: w = 1.0
    # update histogram
    if key in histo: histo[key] += w
    else:            histo[key]  = w    

# calculate free-energy and minimum value
min_fes = 1.0e+15 
for key in histo:
    histo[key] = -KBT_ * math.log(histo[key])
    if(histo[key] < min_fes): min_fes = histo[key]

# total numbers of bins
nbins = 1
for i in range(0, len(nbin)): nbins *= nbin[i]

# print out FES 
log = open(FESFILE_, "w")
# this is needed to add a blank line
xs_old = []
for i in range(0, nbins):
    # get the indexes in the multi-dimensional grid
    key = tuple(get_indexes_from_index(i, nbin))
    # get CV values for that grid point
    xs = get_points(key, gmin, dx)
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
    # print value of CVs
    for x in xs:
        log.write("%12.6lf " % x)
    # print FES
    if key in histo:
       fes = histo[key]-min_fes 
       log.write("   %12.6lf\n" % fes)
    else:
       log.write("       Infinity\n")
log.close()
