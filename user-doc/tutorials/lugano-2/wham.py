import math
import sys

# arguments read from command line
# name of input file
FILENAME_ = sys.argv[1]
# number of BIAS
NBIAS_ = int(sys.argv[2])
# temperature
KBT_ = float(sys.argv[3])
# default parameters for WHAM
# number of WHAM iterations
NWHAM_ = 10000
# convergence thresold
THRES_ = 1.0e-10


def get_wham_weights(nbias, nframes, bias_ts, nwham=NWHAM_, thres=THRES_):
    # find minimum bias
    min_bias = min(bias_ts)
    # initialize weights
    w = []
    for i in range(0, nframes): w.append(1.0)
    # offset and exponential of the bias
    expv = []
    for i in range(0, len(bias_ts)): expv.append(math.exp((-bias_ts[i]+min_bias)/KBT_))
    # initialize Z 
    Z = []
    for j in range(0, nbias): Z.append(1.0)
    # WHAM iterations
    for iii in range(0, nwham):
        # store Z
        Z_old = Z[:]
        # recompute weights
        norm = 0.0
        for i in range(0, len(w)):
            ew = 0.0
            for j in range(0, len(Z)): ew += expv[i*len(Z)+j] / Z[j]
            w[i] = 1.0 / ew
            norm += w[i]
        # normalize weights
        for i in range(0, len(w)): w[i] /= norm
        # recompute Z
        for j in range(0, len(Z)): Z[j] = 0.0
        for i in range(0, len(w)):
            for j in range(0, len(Z)): Z[j] += w[i]*expv[i*len(Z)+j]
        # normalize Z 
        norm = sum(Z)
        for j in range(0, len(Z)): Z[j] /= norm
        # compute change in Z
        eps = 0.0
        for j in range(0, len(Z)): 
            d = math.log(Z[j]/Z_old[j])
            eps += d*d
        # check convergence
        if(eps<thres): break
    # return weights
    return w

# read FILENAME_
bias_ts=[]
for lines in open(FILENAME_, "r").readlines():
    riga=lines.strip().split()
    # skip comment lines
    if(riga[0]=="#!"): continue
    # read bias values
    # umbrella-sampling typical format
    if(len(riga) == NBIAS_+1):
       i0 = 1
       i1 = NBIAS_+1
    # bias exchange typical format
    elif(len(riga) == 2*NBIAS_+1):
       i0 = NBIAS_+1 
       i1 = 2*NBIAS_+1
    # unknown format
    else:
       print(FILENAME_,"format is unknown!")
       exit()
    for i in range(i0, i1):
        bias_ts.append(float(riga[i]))

# number of frames
nframes = len(bias_ts) / NBIAS_

# printout
print("Number of frames::", nframes) 
print("Number of bias entries::", len(bias_ts))

# get wham weights
ws = get_wham_weights(NBIAS_, int(nframes), bias_ts)

# printout WHAM weights to file
print("Weights have been written to weights.dat")
log = open("weights.dat", "w")
for i in range(0, len(ws)): log.write("%32.30lf\n" % ws[i])
log.close()
