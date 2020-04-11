\page Analysis Analysis

PLUMED can be used to analyze trajectories either on the fly during an MD run or via
post processing a trajectory using \ref driver.  A molecular dynamics trajectory is in essence an ordered 
set of configurations of atoms.  Trajectory analysis algorithms are methods that allow us to extract meaningful 
information from this extremely high-dimensionality information.  In extracting this information much of the 
information in the trajectory will be discarded and assumed to be irrelevant to the problem at hand.  For example, 
when we calculate a histogram from a trajectory we throw away all information on the order the frames were visited during the
trajectory.  We instead opt to display a time average that shows the parts of configuration space that were  
visited most frequently.  There are many situations in which this is a reasonable thing to do as we know that
time averages are equivalent to ensemble averages in the long timescale limit and that these average probabilities
of being in different parts of configuration space, \f$P(s)\f$, are thus related to the underlying free
energy, \f$F(s)\f$, via:
\f[
F(s) = - k_B T \ln P(s)
\f]
About the simplest form of analysis 
that PLUMED can perform involves printing information to a file.  PLUMED can output
various different kinds of information to files as described below:

@PRINTANALYSIS@  

The \ref UPDATE_IF action allows you to do more complex things using the above print
commands. As detailed in the documentation for \ref UPDATE_IF when you put any of the above 
actions within an UPDATE_IF block then data will only be output to the file if colvars
are within particular ranges.  In other words, the above printing commands, in tandem 
with \ref UPDATE_IF, allow you to identify the frames in your trajectory that satisfy
some particular criteria and output information on those frames only.

Another useful command is the \ref COMMITTOR command. 
As detailed in the documentation for \ref COMMITTOR this command tells PLUMED (and the underlying 
MD code) to stop the calculation one some criteria is satisfied, alternatively one can use it to keep
track of the number of times a criteria is satisfied.

A number of more complicated forms of analysis can be performed that take a number of frames from 
the trajectory as input.  In all these commands the STRIDE keyword is used to tell PLUMED how 
frequently to collect data from the trajectory.  In all these methods the output from the analysis
is a form of ensemble average.  If you are running with a bias it is thus likely that you may want 
to reweight the trajectory frames in order to remove the effect the bias has on the static behavior
of the system.  The following methods can thus be used to calculate weights for the various trajectory
frames so that the final ensemble average is an average for the canonical ensemble at the appropriate 
temperature.

\section analysisbias Reweighting and Averaging

@REWEIGHTING@

You can then calculate ensemble averages using the following actions.

@GRIDCALC@

For many of the above commands data is accumulated on the grids.  These grids can be further 
analyzed using one of the actions detailed below at some time.  

@GRIDANALYSIS@

As an example the following set of commands instructs PLUMED to calculate the distance between 
atoms 1 and 2 for every fifth frame in the trajectory and to accumulate a histogram from this data
which will be output every 100 steps (i.e. when 20 distances have been added to the histogram).

\plumedfile
x: DISTANCE ATOMS=1,2
h: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 STRIDE=5
DUMPGRID GRID=h FILE=histo STRIDE=100 
\endplumedfile

It is important to note when using commands such as the above the first frame in the trajectory is assumed 
to be the initial configuration that was input to the MD code. It is thus ignored.  Furthermore, if you are 
running with driver and you would like to analyze the whole trajectory (without specifying its length) 
and then print the result you simply call \ref DUMPGRID (or any of the commands above) without a STRIDE 
keyword as shown in the example below. 

\plumedfile
x: DISTANCE ATOMS=1,2
h: HISTOGRAM ARG=x GRID_MIN=0.0 GRID_MAX=3.0 GRID_BIN=100 BANDWIDTH=0.1 STRIDE=5
DUMPGRID GRID=h FILE=histo 
\endplumedfile

Please note that even with this calculation the first frame in the trajectory is ignored when computing the 
histogram.

Notice that all the commands for calculating smooth functions described above calculate some sort of 
average.  There are two ways that you may wish to average the data in your trajectory:

- You might want to calculate block averages in which the first \f$N\f$N frames in your trajectory are
averaged separately to the second block of \f$N\f$ frames.  If this is the case you should use the 
keyword CLEAR in the input to the action that calculates the smooth function.  This keyword is used to 
specify how frequently you are clearing the stored data.

- You might want to calculate an accumulate an average over the whole trajectory and output the average
accumulated at step \f$N\f$, step \f$2N\f$...  This is what PLUMED does by default so you do not need to 
use CLEAR in this case.

\section diag Diagnostic tools

PLUMED has a number of diagnostic tools that can be used to check that new Actions are working correctly: 

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \ref DUMPFORCES </td> <td>Dump the force acting on one of a values in a file.  </td> </tr>
<tr> <td width=5%> \ref DUMPDERIVATIVES </td> <td>Dump the derivatives with respect to the input parameters for one or more objects (generally CVs, functions or biases).</td> </tr>
<tr> <td width=5%> \ref DUMPMASSCHARGE </td> <td>Dump masses and charges on a selected file.</td> </tr>
<tr> <td width=5%> \ref DUMPPROJECTIONS </td> <td>Dump the derivatives with respect to the input parameters for one or more objects (generally CVs, functions or biases).</td> </tr>
</table>

These commands allow you to test that derivatives and forces are calculated correctly
within colvars and functions.  One place where this is very useful is when you are testing whether or
not you have implemented the derivatives of a new collective variables correctly.  So for example if
we wanted to do such a test on the distance CV we would employ an input file something like this:

\plumedfile
d1: DISTANCE ATOMS=1,2
d1n: DISTANCE ATOMS=1,2 NUMERICAL_DERIVATIVES
DUMPDERIVATIVES ARG=d1,d1n FILE=derivatives
\endplumedfile

The first of these two distance commands calculates the analytical derivatives of the distance
while the second calculates these derivatives numerically.  Obviously, if your CV is implemented
correctly these two sets of quantities should be nearly identical.

\section storing Storing data for analysis

All the analysis methods described in previous sections accumulate averages or output diagnostic information on the fly.
That is to say these methods calculate something given the instantaneous positions of the atoms or the instantaneous 
values of a set of collective variables.  Many methods (e.g. dimensionality reduction and clustering) will not work like 
this, however, as information from multiple trajectory frames is required at the point when the analysis is performed.  In other
words the output from these types of analysis cannot be accumulated one frame at time.  When using these methods you must therefore
store trajectory frames for later analysis.  You can do this storing by using the following action:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage COLLECT_FRAMES </td> <td> Collect and store trajectory frames for later analysis with one of the methods detailed below. </td> </tr>
</table>  

\section dissimilaritym Calculating dissimilarity matrices

One of the simplest things that we can do if we have stored a set of trajectory frames using \ref COLLECT_FRAMES is we can calculate the dissimilarity between 
every pair of frames we have stored.  When using the \ref dimred "dimensionality reduction" algorithms described in 
the sections that follow the first step is to calculate this matrix.  Consequently, within PLUMED the following 
command will collect the trajectory data as your simulation progressed and calculate the dissimilarities: 

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage EUCLIDEAN_DISSIMILARITIES </td> <td> Calculate the matrix of dissimilarities between a trajectory of atomic configurations. </td> </tr>
</table>

By exploiting the functionality described in \ref dists you can calculate these dissimilarities in
a wide variety of different ways (e.g. you can use \ref RMSD, or you can use a collection of collective variable
values see \ref TARGET).  If you wish to view this dissimilarity information you can print these quantities 
to a file using:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage PRINT_DISSIMILARITY_MATRIX </td> <td> Print the matrix of dissimilarities between a trajectory of atomic configurations. </td> </tr>
</table>

In addition, if PLUMED does not calculate the dissimilarities you need you can read this information from an 
external file

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage READ_DISSIMILARITY_MATRIX </td> <td> Read a matrix of dissimilarities between a trajectory of atomic configurations from a file. </td> </tr>
</table>
 
N.B. You can only use the two commands above when you are doing post-processing.  

\section landmarks Landmark Selection

Many of the techniques described in the following sections are very computationally expensive to run on large trajectories.
A common strategy is thus to use a landmark selection algorithm to pick a particularly-representative subset of trajectory
frames and to only apply the expensive analysis algorithm on these configurations.  The various landmark selection algorithms
that are available in PLUMED are as follows

@LANDMARKS@

In general most of these landmark selection algorithms must be used in tandem with a \ref dissimilaritym "dissimilarity matrix" object as as follows:

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
data: COLLECT_FRAMES ARG=d1,d2,d3 STRIDE=1
ss1: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=data 
ll2: LANDMARK_SELECT_FPS USE_OUTPUT_DATA_FROM=ss1 NLANDMARKS=300
OUTPUT_ANALYSIS_DATA_TO_COLVAR USE_OUTPUT_DATA_FROM=ll2 FILE=mylandmarks
\endplumedfile

When landmark selection is performed in this way a weight is ascribed to each of the landmark configurations.  This weight is
calculated by summing the weights of all the trajectory frames in each of the landmarks Voronoi polyhedron 
(https://en.wikipedia.org/wiki/Voronoi_diagram).  The weight of each trajectory frame is one unless you are reweighting using the
formula described in the \ref analysisbias to counteract the fact of a simulation bias or an elevated temperature.  If you are reweighting
using these formula the weight of each of the points is equal to the exponential term in the numerator of these expressions.

\section dimred Dimensionality Reduction

Many dimensionality reduction algorithms work in a manner similar to the way we use when we make maps. You start with distances 
between London, Belfast, Paris and Dublin and then you try to arrange points on a piece of paper so that the (suitably transformed) 
distances between the points in your map representing each of those cities are related to the true distances between the cities.  
Stating this more mathematically MDS endeavors to find an <a href="http://en.wikipedia.org/wiki/Isometry">isometry</a> 
between points distributed in a high-dimensional space and a set of points distributed in a low-dimensional plane.  
In other words, if we have \f$M\f$ \f$D\f$-dimensional points, \f$\mathbf{X}\f$, 
and we can calculate dissimilarities between pairs them, \f$D_{ij}\f$, we can, with an MDS calculation, try to create \f$M\f$ projections, 
\f$\mathbf{x}\f$, of the high dimensionality points in a \f$d\f$-dimensional linear space by trying to arrange the projections so that the 
Euclidean distances between pairs of them, \f$d_{ij}\f$, resemble the dissimilarities between the high dimensional points.  In short we minimize:

\f[
\chi^2 = \sum_{i \ne j} w_i w_j \left( F(D_{ij}) - f(d_{ij}) \right)^2
\f]

where \f$F(D_{ij})\f$ is some transformation of the distance between point \f$X^{i}\f$ and point \f$X^{j}\f$ and \f$f(d_{ij})\f$ is some transformation
of the distance between the projection of \f$X^{i}\f$, \f$x^i\f$, and the projection of \f$X^{j}\f$, \f$x^j\f$.  \f$w_i\f$ and \f$w_j\f$ are the weights
of configurations \f$X^i\f$ and \f$^j\f$ respectively.  These weights are calculated using the reweighting and Voronoi polyhedron approaches described in
previous sections.  A tutorial on dimensionality reduction and how it can be used to analyze simulations can be found in the tutorial \ref lugano-5 and in 
the following <a href="https://www.youtube.com/watch?v=ofC2qz0_9_A&feature=youtu.be" > short video.</a>

Within PLUMED running an input to run a dimensionality reduction algorithm can be as simple as:

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
data: COLLECT_FRAMES STRIDE=1 ARG=d1,d2,d3
ss1: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=data 
mds: CLASSICAL_MDS USE_OUTPUT_DATA_FROM=ss1 NLOW_DIM=2
\endplumedfile

Where we have to use the \ref EUCLIDEAN_DISSIMILARITIES action here in order to calculate the matrix of dissimilarities between trajectory frames.
We can even throw some landmark selection into this procedure and perform

\plumedfile
d1: DISTANCE ATOMS=1,2
d2: DISTANCE ATOMS=3,4
d3: DISTANCE ATOMS=5,6
data: COLLECT_FRAMES STRIDE=1 ARG=d1,d2,d3
matrix: EUCLIDEAN_DISSIMILARITIES USE_OUTPUT_DATA_FROM=data
ll2: LANDMARK_SELECT_FPS USE_OUTPUT_DATA_FROM=matrix NLANDMARKS=300
mds: CLASSICAL_MDS USE_OUTPUT_DATA_FROM=ll2 NLOW_DIM=2
osample: PROJECT_ALL_ANALYSIS_DATA USE_OUTPUT_DATA_FROM=matrix PROJECTION=mds
\endplumedfile

Notice here that the final command allows us to calculate the projections of all the non-landmark points that were collected by the action with
label matrix.

Dimensionality can be more complicated, however, because the stress function that calculates \f$\chi^2\f$ has to optimized rather carefully using
a number of different algorithms.  The various algorithms that can be used to optimize this function are described below

@DIMRED@

\section output Outputting the results from analysis algorithms

The following methods are available for printing the result output by the various analysis algorithms:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr> <td width=5%> \subpage OUTPUT_ANALYSIS_DATA_TO_COLVAR </td> <td> Output the results from an analysis using the PLUMED colvar file format. </td> </tr>
<tr> <td width=5%> \subpage OUTPUT_ANALYSIS_DATA_TO_PDB </td> <td> Output the results from an analysis using the PDB file format.</td> </tr>
</table>

Using the above commands to output the data from any form of analysis is important as <b> the STRIDE with which you output the data to a COLVAR or PDB file
controls how frequently the analysis is performed on the collected data </b>.  If you specified no stride on the output lines then PLUMED assumes
you want to perform analysis on the entire trajectory.

If you use the above commands to output data from one of the \ref landmarks algorithms then only the second will give you information on the 
atomic positions in your landmark configurations and their associated weights.  The first of these commands will give the values of the colvars
in the landmark configurations only.  If you use the above commands to output data from one of the \ref dimred algorithms then 
\ref OUTPUT_ANALYSIS_DATA_TO_COLVAR will give you an output file that contains the projection for each of your input points.  \ref OUTPUT_ANALYSIS_DATA_TO_PDB
will give you a PDB that contains the position of the input point, the projections and the weight of the configuration.

A nice feature of plumed is that when you use \ref landmarks algorithms or \ref dimred algorithms the output information is just a vector of 
variables.  As such you can use \ref HISTOGRAM to construct a histogram of the information generated by these algorithms.

