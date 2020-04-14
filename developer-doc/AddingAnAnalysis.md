\page AddingAnAnalysis Implementing analysis methods in PLUMED

Information on implementing methods that perform analysis on stored trajectory information eg dimensionality reduction

Implementing methods for analysing trajectory data is more complex than implementing collective variables.
Consequently it is difficult to write a step by step guide like those we have written on implementing \ref AddingAColvar "colvars" or
\ref AddingAFunction "functions".  Hence, this document tries to explain the things that we have considered and the way these 
have been incorporated into the PLMD::analysis::AnalysisBase and PLMD::analysis::AnalysisWithDataCollection abstract base classes.  
Hopefully this will provide some insight into our rationale in writing these parts of the code and will help you to understand how 
any new analysis method can be implemented in the PLUMED code in a way that exploits those features that are already there.

\section overview An overview of analysis methods in PLUMED

There are two distinct ways in which one may wish to perform some form of analysis on a molecular dynamics trajectory.  In the first method some quantity
is calculated for each of the atoms in the trajectory in each of the frames and this is then averaged over the whole trajectory.  Velocity 
autocorrelation functions or mean squared displacements are examples of forms of analysis of this type.   The methods implemented in PLMD::analysis are 
not of this type.  These methods are designed to collect set of snapshots of the trajectory and to perform some analysis of these snapshots.
These trajectory snapshots might be the values of a particular set of collective variables for each of the frames, they might be the 
instantaneous positions of the atoms in each frame or they might be some combination of the above.  The assumption then when running one of these
analysis methods is that a representation (or snapshot) will be collected intermittently from the trajectory and then once a sufficiently large
 collection of these snapshots are collected they will be analysed. 

\section antraj Analysis on the fly

It is important to remember that PLUMED is primarily a code for performing biased molecular dynamics.  The code's original purpose was
to be a plugin that could be easily added to a number of different molecular dynamics engines that allowed additional forces to be 
incorporated when integrating the equations of motion.  The command line tools that allows one to analyse trajectories during post
processing were added later as an afterthought.  This consideration is particularly important when considering analysis algorithms 
in this code because most analysis codes that are used in the community read in the trajectory and to do all their calculations 
during post processing.  The analysis functions that have been implemented in PLUMED can be used to post-process trajectories - you
simply make use of the command line tool driver - however, they can also be used to analyse trajectories on the fly.  We believe this 
is useful for a number of reasons:

- Computers are becoming more powerful so it is possible to run simulations for much longer.  At the same time, however, hard drives
space is at a premium and it is thus difficult to store these large trajectories for post-processing.  If trajectories can be analysed
on the fly this presents less of a problem.
- A number of free energy methods (e.g. Gaussian mixture umbrella sampling, adaptive umbrella sampling and reconnaissance metadynamics)
work by performing a sophistacated analysis of the trajectory and then using the result from this analysis to design a bias for the 
dynamics. 
- Analysis methods implemented in PLUMED can take advantage of the many different collective variables that we have implemented in 
this code and are thus extremely flexible implementations of these techniques.
- Analysis methods implemented in PLUMED make use of the PLUMED input syntax, which hopefully allows users already familiar with 
PLUMED to get to grips with using these tools more rapidly. 

\section genphil General Phillosopy

The aim in the PLMD::analysis and PLMD::dimred modules is to write code that is very flexible and that allows the user a great deal of 
flexibility in the input.  For example we would like to be able to write the code so that a user can collect data from the trajectory and 
then at analysis time they can:

- Select a subset of landmark points from the stored data
- Generate projections of these landmark points using sketch-map
- Project the remaining non-landmark points using the sketch-map projection generated and construct a histogram as a function of the sketch-map coordinates.

Furthermore, we would like to be able to do all the above with a minimum of syntax for the user and with a minimum amount of code to maintain.
This is why the analysis class is structured in the way it is.  The general idea is that one PLMD::analysis::AnalysisWithDataCollection object 
collects the trajectory data as the simulation progresses (for the example above this would be an object of type PLMD::analysis::EuclideanDissimilarityMatrix).
Then when it is time to analyse a chain of analysis objects are run on the data collected by the PLMD::analysis::AnalysisWithDataCollection.
There are thus two types of analysis actions:

- Those that inherit from PLMD::analysis::AnalysisWithDataCollection - these can collect data from a trajectory
- Those that inherit from PLMD::analysis::AnalysisBase - these cannot collect data from a trajectory.  They get their input from another PLMD::analysis::AnalysisBase Action

If you look at the code most of the routines in the PLMD::analysis and PLMD::dimred modules do not inherit from PLMD::analysis::AnalysisWithDataCollection.
In fact the only ones where this is an option that users really see in the manual are PLMD::analysis::Histogram and PLMD::analysis::EuclideanDissimilarityMatrix
(many of the other analysis actions that inherit from PLMD::analysis::AnalysisWithDataCollection only do so because this allows us to write straightforward 
regression tests for this part of the code).  The remaining analysis actions inherit from PLMD::analysis::AnalysisBase because in general they require some 
dissimilarity information in order to function and because this information is calculated within PLMD::analysis::EuclideanDissimilarityMatrix.  
Consequently, if you are writing a new analysis class it should probably not inherit from PLMD::analysis::AnalysisWithDataCollection.

\section storing Storing trajectory data for analysis

As discussed in the previous section the methods in PLMD::analysis all work by storing a snapshot of the trajectory every \f$N\f$ steps
and by then running an analysis method every \f$M\f$ steps, where \f$M\f$ is some substantial multiple of \f$N\f$.  This intermittent
storing of trajectory frames and occasional analysis of the trajectory is all looked after within the PLMD::analysis::AnalysisWithDataCollection
abstract base class.  Users can then set of a chain of analysis actions on the stored data by using actions that inherit from
PLMD::analysis::AnalysisBase.   Any class inheriting from PLMD::analysis::AnalysisBase must have a method within it named performAnalysis(),
which will actually perform the analysis on the stored trajectory frames.  When implementing a new analysis method the majority of your
development time will probably be spent implementing some part of this performAnalysis method.

\section reweight Reweighting trajectories

As discussed in previous sections PLUMED is primarily a code for doing biased molecular dynamics simulations.  This bias is used to force rare events
to occur in the short time scales that are accessible within a molecular dynamics simulation.  Oftentimes when analysing a trajectory we would like
to remove the contribution of the bias and to reweight the simulation so as to get the underlying free energy surface.  When performing any analysis of
the trajectory one may similarly wish to remove the effect of the bias and to consider what weight each sampled point would have had if it had been sampled
in accordance with the underlying canonical probability distribution function.  This process of reweighting points - of ascribing a weight to each snapshot
that discounts the effect of the simulation bias - is again looked after within PLMD::analysis::AnalysisWithDataCollection.  If you wish to take these
weights into account in your analysis method you should use the method PLMD::analysis::AnalysisBase::getWeight to access them.  Obviously, if you have
no simulation bias on the system then each point will have a weight of one and this will be the weight returned by PLMD::analysis::AnalysisBase::getWeight.

\section output Outputting data files

The fact that PLMD::analysis::AnalysisWithDataCollection can be used to run trajectory analysis in either post processing or on the fly during a trajectory means
that this class must also look after a number of things.  For example one might wish to perform multiple analyses of the trajectory 
during a simulation.  Obviously, you do not want to overwrite the output file from your first analysis when you perform the second 
analysis of the trajectory.  In addition, you do not want to overwrite files from earlier runs if you choose to rerun your analysis 
in a directory where you had already run an earlier calculation.  For these reasons whenever you wish to read in the name of an output file 
you should use the following code to make sure that any old files are backed up on restart:

\verbatim
if( !getRestart() ){ OFile ofile; ofile.link(*this); ofile.setBackupString("analysis"); ofile.backupAllFiles(fname); }
\endverbatim 

where fname is the name of your output file. On top of this when you open an output file in your analysis method you should use the following 
set of commands:

\verbatim
OFile gfile; gfile.link(*this);
gfile.setBackupString("analysis");
gfile.open( ofilename.c_str() ); 
\endverbatim

The second line ensures that files are named analysis.0.ofilename, analysis.1.ofilename and so on. Having said that it is probably best 
to not write routines to output data in analysis classes and to simply ensure that you can pass any output data from your method to the 
PLMD::analysis::OutputColvarFile and PLMD::analysis::OutputPDBFile methods.  If you have done everything properly these classes should be
able to interact with your new class through methods that are in PLMD::analysis::AnalysisBase.

\section metrics Calculating dissimilarity matrices

One of the most basic ways in which you can analyse a collection of trajectory snapshots is to calculate the matrix of dissimilarities between each of the pairs
of trajectories frames.  In plumed this is done in by the class PLMD::analysis::EuclideanDissimilarityMatrix.  Notice as well that this class makes full use
of the PLMD::reference module that is discussed in \ref AddingAMetric and so this class alone allows you to calculate the dissimilairity between any pair of 
trajectory frames in a wide variety of different ways.  In addition, you can use PLMD::analysis::ReadDissimilarityMatrix to get the dissimilarities from an 
input file.  There should thus be little need to implement new objects for calculating dissimilarity matrices -- if you really feel you need something other than
what is already there or that is not implementable by \ref AddingAMetric then you are doing a calculation that is very unusual.

\section landmarks Landmark selection algorithms

Many analyses methods scale poorly with the number of trajectory frames that you wish to analyse.  This happens in part because you need to compute the matrix
of pairwise disimiarities (a square matrix in which the number of columns is equal to the number of trajectory frames) but also because you then have to 
do some algebra involving this matrix.  To alleviate these problems a common strategy is to perform the analysis on a set of so-called landmark frames and to 
then project the non-landmark snapshots from the trajectory using some out-of-sample extension of your analysis algorithm.  Classes that inherit from
PLMD::analysis::LandmarkSelectionBase are implementations of the various landmark selection algorithms that are commonly employed.  If your favourite landmark
selection algorithm is not there you may choose to implement a new landmark selection algorithm by writing a new PLMD::Action that inherits from this class.
Within this new class you need to only define a registerKeywords function, a constructor and a method that actually selects the landmarks that will be a function
that must be called selectLandmarks.  Once you are satisfied that you want frame \f$k\f$ in your landmark set then you can select this using 
PLMD::analysis::LandmarkSelectionBase::selectFrame.  Please note that you can get the number of landmarks to select by calling:

\verbatim
getNumberOfDataPoints()
\endverbatim

If you would like to get the total number of frames from which you can get your subset of landmarks you should call:

\verbatim
mydata->getNumberOfDataPoints()
\endverbatim

Lastly, if you are using a method, which like PLMD::analysis::FarthestPointSampling uses the distances between your input points in some way, you should 
add something akin to:

\verbatim
if( !dissimilaritiesWereSet() ) error("dissimilarities have not been calcualted in input action");
\endverbatim

in your constructor as this will ensure that users are told what they are doing wrong if they do not specify how to calculate distances between points in the
input.  To then get the dissimilirity between input point \f$i\f$ and input point \f$j\f$ use:

\verbatim
mydata->getDissimilarity( landmarks[i], k );
\endverbatim

Calling PLMD::analysis::AnalysisBase::getDissimilarity will give you the distance between a pair of landmarks, which is not what you need.

\section dimred Dimensionality reduction

The aim when writing any large code such as PLUMED is to minise the number of lines of code as fewer lines of code means fewer bugs on average.
Hence, as explained in other sections of this developer manual, all the object oriented programming, inheritance and polymorphism.  Given this
consider how we would go about implementing a library of dimensionality reduction algorithms.  In LLE, ISOMAP, sketch-map or MDS the aim is to
generate a low-dimensional projection of some set of high-dimensional data points.  For all these methods we can use the same code to to store 
the high and low dimensional points and to output this data to a file.  In fact the only things that differ in these various different methods are
the ways in which the dissimilarities between the high-dimensional points are calculated and the manner in which the low-dimensional projections
are generated.  We have already discussed how PLUMED calculates matrices of dissimilarities between points using PLMD::analysis::EuclideanDissimilarityMatrix
and how one can think about introduce new methods of calculating dissimilarities.  Priting to files meanwhile can be looked after by PLMD::analysis::OutputColvarFile 
and PLMD::analysis::OutputPDBFile.  Furthermore, if dimensionality reduction classes are written properly it is even possible to pass the projections
generated by them to PLMD::analysis::Histogram and to output histograms and free energy surfaces as a function of the low-dimensional coordinates.
As such in PLMD::dimred::DimensionalityReductionBase and its daughter classes we are thus solely concerned with calculating projections of data points.  
The dissimilarities between the input high dimensional frames will always be calculated by some method akin to PLMD::analysis::EuclideanDissimilarityMatrix.  
PLMD::dimred::DimensionalityReductionBase inherits from PLMD::analysis::AnalysisBase and expects another PLMD::analysis::AnalysisBase object as input.  
Furthermore, with the exception of PLMD::dimred::ClassicalMultiDimensionalScaling the input PLMD::AnalysisBase must be an 
PLMD::dimred::DimensionalityReductionBase as initial guesses must be suggested when using an iterative optimization algorithm such as PLMD::dimred::SMACOF.

Much of the work of dimensionality reduction done in the base PLMD::dimred::DimensionalityReductionBase class.  When implementing any new
dimensionality reduction algorithm your focus will be on writing a PLMD::dimred::calculateProjections routines.  This function will take as
input the matrix of pairwise dissimilarities between points (be they landmark points or otherwise) and is expected to return a matrix 
containing the projections of the high-dimensional data points.  Depending on the stress function you minimise to find projections you 
may also have to implement PLMD::dimred::DimensionalityReductionBase::setTargetDistance and PLMD::dimred::DimensionalityReductionBase::calculateStress
functions.  This is necessary with PLMD::dimred::SketchMapBase for example because sketch-map uses transformed versions of the dissimilarities
and distances in its stress function.  These two routines are used within PLMD::dimred::ProjectNonLandmarkPoints which will product optimal projections of 
points that were not selected as landmarks.   

\section cluster Clustering trajectory data

There are currently no clustering methods implemented in the PLMD::analysis module.  This section is thus here to explain how I (Gareth Tribello) imagined one
might go about implementing these methods.  Again there are many commonalities between methods such as kmeans, Gaussian mixture models and so on, which should 
be thought about when constructing an abstract base class.  Furthermore, this abstract base class should (like PLMD::analysis::DimensionalityReductionBase) 
inherit from PLMD::analysis::AnalysisBase and be implemented in a way that allows one to exploit the use of landmark selection algorithms that inherit from 
PLMD::analysis::LandmarkSelectionBase and the measure actions such as PLMD::analysis::EuclideanDissimilarityMatrix.  It should also be able to work with the 
weights you get by reweiting the trajectories with the bias and so on.  I don't think that you should have this inheriting from 
PLMD::analysis::AnalysisWithDataCollection as I believe you want the clustering to work with projections of data generated by dimensionality reduction algorithms.
Obviously, if you are thinking of adding methods to cluster trajectory frames within PLUMED please feel free to get in touch with me (gareth.tribello\@gmail.com).  
I will be more than happy to discuss these ideas with you.

