import glob
import numpy as np
# Here are some numbers you will need to change if you run this script on grids generated in different contexts
nquantities = 1          # Number of quanities that have been evaluated on the grid
grid_dimension = 1       # Number of collective variables that you provided using the ARG keyword
filename = "myhist.dat"  # The name you specified the data to output to in the DUMPGRID command 
# Function to read in histogram data and normalization
def readhistogram( fname ) :
        # Read in the histogram data
        data = np.loadtxt( fname )
        with open( filename, "r" ) as myfile :
                for line in myfile :
                        if line.startswith("#! SET normalisation") : norm = line.split()[3]
        return float(norm), data
# Read in the grid file header to work out what fields we have
with open( filename, "r" ) as myfile :
        for line in myfile :
                if line.startswith("#! FIELDS") : fieldnames = line.split()
# Check if derivatives have been output in the grid by investigating the header
nextg = 1
if len(fieldnames)>(2+grid_dimension+nquantities) :
        nextg = 1 + grid_dimension
        assert len(fieldnames)==(2+grid_dimension + nquantities*nextg)
# Read in a grid
norm, griddata = readhistogram( filename )
norm2 = norm*norm
# Create two np array that will be used to accumulate the average grid and the average grid squared
average = np.zeros((nquantities, len(griddata[:,0])))
average_sq = np.zeros((nquantities, len(griddata[:,0])))
for i in range(0,nquantities) :
        average[i,:] = norm*griddata[:,nquantities+i*nextg]
        average_sq[i,:] = norm*griddata[:,nquantities+i*nextg]*griddata[:,nquantities+i*nextg]
# Now sum the grids from all all the analysis files you have
for filen in glob.glob( "analysis.*." + filename ) :
        tnorm, newgrid = readhistogram( filen )
        norm = norm + tnorm
        norm2 = norm2 + tnorm*tnorm
        for i in range(0,nquantities) :
                average[i,:] = average[i,:] + tnorm*newgrid[:,nquantities+i*nextg]
                average_sq[i,:] = average_sq[i,:] + tnorm*newgrid[:,nquantities+i*nextg]*newgrid[:,nquantities+i*nextg]
# Compte the final average grid
average = average / norm
# Compute the sample variance for all grid points
variance = (average_sq / norm) - average*average
# Now multiply by bessel correction to unbias the sample variance and get the population variance
variance = ( norm /(norm-(norm2/norm)) ) * variance
# And lastly divide by number of grids and square root to get an error bar for each grid point
ngrid = 1 + len( glob.glob( "analysis.*." + filename ) )
errors = np.sqrt( variance / ngrid )
# Calcualte average error over grid and output in header
for i in range(0,nquantities) :
        mean_error = sum(errors[i,:]) / len(errors[i,:])
        print("# Average error for " + str(i+1) + "th averaged function on grid equals ", mean_error )
# Output the final average grid
for i in range(0,len(griddata[:,0])) :
        for j in range(0,grid_dimension) : print( griddata[i,j], end=" " )
        for j in range(0,nquantities) : print( average[j,i], errors[j,i], end=" " )
        print()
