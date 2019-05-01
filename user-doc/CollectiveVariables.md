\page colvarintro Collective Variables

Chemical systems contain an enormous number atoms, which, in most cases makes it simply impossible for
us to understand anything by monitoring the atom positions directly.  Consequently,
we introduce Collective variables (CVs) that describe the chemical processes we are
interested in and monitor these simpler quantities instead.  These CVs are used in many of the methods
implemented in PLUMED - there values can be monitored using \ref PRINT, \ref Function of them can be calculated
or they can be analyzed or biased using the \ref Analysis and \ref Bias "Biasing" methods implemented in PLUMED.
Before doing any of these things however we first have to tell PLUMED how to calculate them.

The simplest collective variables that are implemented in PLUMED take in a
set of atomic positions and output one or multiple scalar CV values.  Information on these variables is given on the page entitled 
\ref Colvar while information as to how sets of atoms can be selected
can be found in the pages on \ref Group.  Please be aware that PLUMED contains implementations of many other collective variables 
but that the input for these variables may be less transparent when it is first encountered.
In particular, the page on \ref dists describes the various ways that you can calculate the distance from a particular reference
configuration.  So you will find instructions on how to calculate the RMSD distance from the folded state of a protein here.
Meanwhile, the page on \ref Function describes the various functions of collective variables that can be used in the
code.  This is a very powerful feature of PLUMED as you can use the \ref Function commands to calculate any function or 
combination of the simple collective variables listed on the page \ref Colvar.  Lastly the page on \ref mcolv describes MultiColvars.  
MultiColvars allow you to use many different colvars and allow us to
implement all these collective variables without a large amount of code.  For some things (e.g.
\ref DISTANCES GROUPA=1 GROUPB=2-100 LESS_THAN={RATIONAL R_0=3}) there are more computationally efficient options available in plumed
(e.g. \ref COORDINATION).  However, MultiColvars are worth investigating as they provide a flexible syntax for many quite-complex CVs.

- \subpage Group
- \subpage Colvar
- \subpage dists
- \subpage Function
- \subpage mcolv
- \subpage contactmatrix

\page Colvar CV Documentation

The following list contains descriptions of a number of the colvars that are currently implemented in PLUMED.

@COLVAR@

\page dists Distances from reference configurations

One colvar that has been shown to be very successful in studying protein folding is the distance between the instantaneous configuration
and a reference configuration - often the structure of the folded state.  When the free energy of a protein is shown as a function
of this collective variable there is a minima for low values of the CV, which is due to the folded state of the protein.  There is 
then a second minima at higher values of the CV, which is the minima corresponding to the unfolded state.

A slight problem with this sort of collective variable is that there are many different ways of calculating the distance from a 
particular reference structure.  The simplest - adding together the distances by which each of the atoms has been translated in
going from the reference configuration to the instantaneous configuration - is not particularly sensible.  A distance calculated
in this way does not neglect translation of the center of mass of the molecule and rotation of the frame of reference.  A common practice
is thus to remove these components by calculating the \ref RMSD distance between the reference and instantaneous configurations.
This is not the only way to calculate the distance, however.  One could also calculate the total amount by which a large number 
of collective variables change in moving from the reference to the instantaneous configurations.  One could even combine RMSD distances
with the amount the collective variables change.  A full list of the ways distances can be measured in PLUMED is given below:

@DCOLVAR@

These options for calculating distances are re-used in a number of places in the code.  For instance they are used in some of the 
analysis algorithms that are implemented in PLUMED and in \ref PATH collective variables. 
Notice that most of these actions read the reference configuration from a PDB file. Be sure
you understand how to format properly a PDB file to use used in PLUMED (see \ref pdbreader).

\page mcolv MultiColvar 

Oftentimes, when you do not need one of the collective variables described elsewhere in the manual, what you want instead is a 
function of a distribution of collective variables of a particular type.  In other words, you would like to calculate a
function something like this:
\f[
s = \sum_i g[f(\{X\}_i)]
\f]
In this expression \f$g\f$ is a function that takes in one argument and \f$f\f$ is a function that takes a set of atomic positions
as argument. The symbol \f$\{X\}_i\f$ is used to indicate the fact that the function \f$f\f$ is evaluated for a number of different
sets of atoms.  If you would just like to output the values of all the various \f$f\f$ functions you should use the command \ref DUMPMULTICOLVAR

This functionality is useful if you need to calculate a minimum distance or the number of coordination numbers greater than a 3.0.  
To avoid duplicating the code to calculate an angle or distance many times and to make it easier to implement very complex collective 
variables PLUMED provides these sort of collective variables using so-called MultiColvars.  MultiColvars are named in this way because a single
PLUMED action can be used to calculate a number of different collective variables.  For instance the \ref DISTANCES
action can be used to calculate the minimum distance, the number of distances less than a certain value, the number of
distances within a certain range... A more detailed introduction to multicolvars is provided in this 
<a href="http://www.youtube.com/watch?v=iDvZmbWE5ps">10-minute video</a>. Descriptions of the various multicolvars
that are implemented in PLUMED 2 are given below: 

@MCOLVAR@  

To instruct PLUMED to calculate a multicolvar you give an instruction that looks something like this:

\verbatim
NAME <atoms involved> <parameters> <what am I calculating> TOL=0.001 LABEL=label
\endverbatim

Oftentimes the simplest way to specify the atoms involved is to use multiple instances of the ATOMS keyword 
i.e. ATOMS1, ATOMS2, ATOMS3,...  Separate instances of the quantity specified by NAME are then calculated for 
each of the sets of atoms.  For example if the command issued contains the following:

\plumedfile
DISTANCES ATOMS1=1,2 ATOMS2=3,4 ATOMS3=5,6
\endplumedfile

The distances between atoms 1 and 2, atoms 3 and 4, and atoms 5 and 6 are calculated. Obviously, generating 
this sort of input is rather tedious so short cuts are also available many of the collective variables. 
These are described on the manual pages for the actions.
 
After specifying the atoms involved you sometimes need to specify some parameters that required in the 
calculation.  For instance, for \ref COORDINATIONNUMBER - the number of atoms in the first coordination
sphere of each of the atoms in the system - you need to specify the parameters for a \ref switchingfunction
that will tell us whether or not an atom is in the first coordination sphere.  Details as to how to do this
are provided on the manual pages.  

One of the most important keywords for multicolvars is the TOL keyword.  This specifies that terms in sums 
that contribute less than a certain value can be ignored.  In addition, it is assumed that the derivative
with respect to these terms are essentially zero.  By increasing the TOL parameter you can increase the speed 
of the calculation.  Be aware, however, that this increase in speed is only possible because you are lowering the 
accuracy with which you are computing the quantity of interest.

Once you have specified the base quantities that are to be calculated from the atoms involved and any parameters
you need to specify what function of these base quantities is to be calculated.  For most multicolvars you can calculate 
the minimum, the number less than a target value, the number within a certain range, the number more than a target
value and the average value directly.  

\section multicolvarfunction MultiColvar functions

It is possible to use multicolvars to calculate complicated collective variables by exploiting the fact that the output
from one multicolvar can be used as input to a second multicolvar.  One simple way of exploiting this functionality is to
filter the atoms based on the value they have for a symmetry function.  For example you might want to consider only those
atoms that with a \ref COORDINATIONNUMBER higher that a certain threshold when calculating some particularly expensive symmetry
function such as \ref Q6.  The following methods can thus all be used to filter the values of multicolvars in this way:

@MFILTERS@

An alternative way of filtering atoms is to consider only those atoms in a particular part of the simulation box.  This can be
done by exploiting the following methods

@VOLUMES@

The idea with these methods is that function of the form:
\f[
s = \sum_i w(\{X\}_i) g[f(\{X\}_i)]
\f]
can be evaluated where once again \f$g\f$ is a function with one argument and \f$g\f$ is a function of a set of atomic positions.  
The difference from the more general function described earlier is that we now have a weight \f$w\f$ which is again a function of the
atomic positions.  This weight varies between zero and one and it is this weight that is calculated in the list of filtering methods
and volume methods described in the lists above.  

In addition to these volume and filtering methods it is also possible to calculate the local average of a quantities in the manner 
described in \cite dellago-q6 using the \ref LOCAL_AVERAGE method.  Furthermore, in many cases \ref Q6, \ref MOLECULES and 
\ref PLANES the symmetry function being evaluated is a vector.  You can thus construct a variety of novel collective variables by
taking dot products of vectors on adjacent atoms as described below: 

@MCOLVARF@ 

The final set of functions that you can apply on multicolvars are functions that transform all the colvars calculated using a 
multicolvar using a function.  This can be useful if you are calculating some complicated derived quantity of some simpler 
quantity.  It is also useful if you are calculating a Willard Chandler surface or a histogram.  The actions that you can use to 
perform these transforms are:

@MTRANSFORMS@

\section multicolvarbias MultiColvar bias

There may be occasions when you want add restraints on many collective variables. For instance if you are studying a cluster
you might want to add a wall on the distances between each of the atoms and the center of mass of the cluster in order to
prevent the cluster subliming.  Alternatively, you may wish to insist that a particular set of atoms in your system all have a 
coordination number greater than 2.  You can add these sorts of restraints by employing the following biases, which all act 
on the set of collective variable values calculated by a multicolvar.  So for example the following set of commands:

\plumedfile
COM ATOMS=1-20 LABEL=c1
DISTANCES GROUPA=c1 GROUPB=1-20 LABEL=d1
UWALLS DATA=d1 AT=2.5 KAPPA=0.2 LABEL=sr
\endplumedfile

creates the aforementioned set of restraints on the distances between the 20 atoms in a cluster and the center of mass of the cluster.

The list of biases of this type are as follows:

@MCOLVARB@

Notice that (in theory) you could also use this functionality to add additional terms to your force field or to implement your 
force field.

\section usingbase Extracting all the base quantities

There may be occasions where you want to get information on all the individual colvar values that you have calculated.
For example you might want to output the values of all the coordination numbers calculated by a \ref COORDINATIONNUMBER 
action.  You can thus use the following command to extract this sort of information, \ref DUMPMULTICOLVAR.

\page contactmatrix Exploiting contact matrices

A contact matrix is an \f$N \times N\f$ matrix in which the \f$i\f$th, \f$j\f$th element tells you whether or not the \f$i\f$th
and \f$j\f$th atoms/molecules from a set of \f$N\f$ atoms/molecules are adjacent or not.  There are various ways of defining
whether a pair of atoms/molecules are adjacent or not.  For example we can say two atoms are adjacent if the distance between
them is less than some cutoff.  Alternatively, if we have a have a pair of molecules, we might state they are adjacent if their
centers of mass are within a certain cutoff and if the two molecules have the same orientation.  Two electronegative atoms
might be said to be adjacent if there is a hydrogen bond between them.  For these reasons then PLUMED contains all of the 
following methods for calculating an adjacency matrix 

@MATRIX@

Once you have calculated an adjacency matrix you can then perform any one of the following operations on this object in order
to reduce it to a scalar number or a set of connected components.

@MATRIXF@

If the function you have chosen reduces your contact matrix to a set of connected components you then need a method to convert 
these connected components into a scalar number or to output this information to a file.  The various things that you can do
with a set of connected components are listed below:

@CONCOMP@

