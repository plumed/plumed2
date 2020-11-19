\page AddingAMetric Implementing methods for calculating the distances between pairs of configurations 

To implement a new method for calculating the distance between a pair of trajectory frames you will need to work with the
PLMD::reference module.  This module is used in many parts of PLUMED (e.g. path collective variables, field-cvs and analysis methods).
Consequently, if you implement your distance measure using the functionality in this class you will be able to use it in a 
wide variety of difffernt contexts.  As always if you would like us to incorporate your measure in the release version of 
PLUMED you will need to write at least one regression test for it.  Details on how to write regression tests are provided 
here: \ref regtests  

The classes in PLMD::reference allow one to calculate the distance between pairs of trajectory frames.  The reason that this is 
a module rather than a single method is that there are a wide variety of different ways of calculating the distance between
two frames.  One could, for example, calculate the RMSD distance between the two frames.  Alternatively, one can calculate 
a large number of collective variables and compare the values these variables have in the two configurations.  Lastly, one
could somehow combine some element of RMSD calculation with the calculation of some set of collective variables.  As with
so many things the way in which one wishes to calculate the distance between two configurations will depend on the problem
one is endeavoring to solve.  The aim of the PLMD::reference module is thus to provide the user unlimited flexibility in the
way that the matrix of distances can be calculated in analysis methods such as sketch-map or in biasing methods such as path
cvs.  I say unlimited because, although the particular distance measure the user needs many not currently be implemented in 
PLUMED, he/she always has the option to implement this new distance metric in the reference module and then access it in the
full range of analysis and biasing tools that are already available in the code.  The following provides instructions for 
implementing a new way of calculating the distance between a pair of trajectory frames in the code using the PLMD::reference 
module.  As always once your new method is implemented you can access it in all the places where pairwise distances are employed
because of the way that PLUMED exploits inheritance and polymorphism.   

\section adding Creating your measure

To create a new way of measuring the distances between pairs of atoms you must write a new class in the reference directory.
An example declaration for such a class is given below:

\verbatim
class OptimalRMSD : public RMSDBase {
private:
  bool fast;
  RMSD myrmsd;
public:
  OptimalRMSD(const ReferenceConfigurationOptions& ro);
  void read( const PDB& );
  double calc( const std::vector<Vector>& pos, const bool& squared );
};
\endverbatim

In this case we are inheriting from PLMD::RMSDBase but the particular class you will want to inherit from will depend on the
manner in which the distance is calculated.  To then ensure that your new measure can be used throughout the code you need
to include the MetricRegister.h file and the following line:

\verbatim
PLUMED_REGISTER_METRIC(OptimalRMSD,"OPTIMAL")
\endverbatim

Once again dynamic polymorphism is exploited again here.  With this example the command above ensures that PLMD::OptimalRMSD objects
are used whenever the command METRIC=OPTIMAL is found in the PLUMED input.

Your new class must contain a constructor, a method to read in the configuration from a pdb file (read) and a method to calculate
the distance between the reference configuration and the input position (calc).  Please be aware that the number of arguments
to this method will change depending on which base class you inherit from in creating your new measure object.

The inheritance structure of these routines is rather complicated looking.  In essence, however, the base class underlying all
these classes in PLMD::ReferenceConfiguration.  There are then two classes PLMD::ReferenceArguments and PLMD::ReferenceAtoms that
can be multiply inherited in any derived classes.  PLMD::ReferenceArguments provides tools for dealing with distance measures that
involve colvars. PLMD::ReferenceAtoms provides tools for dealing with distance measures that involve the positions of atoms.
Base classes such as PLMD::SingleDomainRMSD and PLMD::RMSDBase are there essentially so that we can use particular sets of metrics
in secondary structure variables and the RMSD class respectively.  Essentially within these two objects the RMSD is calculated by
calling the calc member of the abstract base class PLMD::SingleDomainRMSD and PLMD::RMSDBase.  This allows one to use the
appropriate set of measures within these particular functions.

\section args Dealing with colvars

There are essentially three ways in which you might wish to calculate the distance between two sets of reference colvar values:

- You will calculate the euclidean distance using pythagoras theorem
- You calculate the normalised euclidean distance in which pythagoras theorem is again used but each pair of components is given a separate weight
- You calculate the Mahalonobis distance in which the distance is calculated as \f$x^T M x\f$ where \f$x\f$ is the displacement vector and \f$M\f$ is a matrix.

These three methods are all implemented within PLUMED and are in the classes PLMD::EuclideanDistance, PLMD::NormalizedEuclideanDistance and PLMD::MahalanobisDistance
respectively.  If you look in these three classes you will note that there is a very small amount of code in each of them.  Essentially the calculation of the
distance and the reading in of the PDB file are looked after by the methods PLMD::ReferenceArguments::calculateArgumentDistance and
PLMD::ReferenceArguments::readArgumentsFromPDB respectively.  To reuse these functionalities you need to add the command:

- hasmetric=true in the constructor if you want to use the Mahalonobis distance in your new metric
- hasweights=true in the constructor if you want to use the normalised euclidean distance in your new metric

If you want to use the euclidean distance you can use PLMD::ReferenceArguments::calculateArgumentDistance and PLMD::ReferenceArguments::readArgumentsFromPDB
without any further instructions.  Notice these methods will still work if you use a combination of atom positions and colvars in your measure.  They will
give the part of the measure due to the arguments - you will be required to do some further calculation to get the bit involving the atoms.

\section atoms Dealing with atoms

All the distance measures that we work with involve the positions of the atoms in the two configurations.  If we work with colvars we
just do the calculation of the distance from the positions of the atoms in an indirect way - we calculate some intermediary quantities
from the positions of the atoms and then calculate the set of difference between these intermediary quantities.  There are cases,
however, such as when we are working with RMSD distances where it is useful to calculate the distance directly from the set of atomic
positions.  That is to say there are cases where these intermediary quantities are not useful.  If you are implementing such a measure
you will need to write a class that inherits from PLMD::ReferenceAtoms either directly or indirectly.

Within the read method you can read the atoms in the PDB file by using the method PLMD::ReferenceAtoms::readAtomsFromPDB.  Within calc you can then
access the positions of these read in atoms by using the method PLMD::ReferenceAtoms::getReferencePositions() or by using PLMD::ReferenceAtoms::getReferencePosition.
It is useful to think carefully about how you can reuse other parts of the code when implementing new reference methods.  As an example notice
how PLMD::OptimalRMSD and PLMD::SimpleRMSD make extensive use of the PLMD::RMSD class.  Similarly notice how the PLMD::OptimalRMSD,
PLMD::SimpleRMSD and PLMD::DRMSD classes are reused in PLMD::MultiDomainRMSD.

\section AddingAMeasureDocs Adding documentation for your measure

To test whether you have implemented you new measure correctly you should implement a small Action that calculates the
distance between the instantaneous positions of the atoms and some read in reference configuration.  The methods
PLMD::colvar::RMSD, PLMD::colvar::DRMSD, PLMD::colvar::MultiRMSD and PLMD::function::Target perform these functions for
the set of measures that are currently implemented within the reference module.  You will need to do something similar
in your test action.

You will notice that the documentation in these files starts with the line:

\verbatim
//+PLUMEDOC DCOLVAR TARGET
\endverbatim

The DCOLVAR tag is important here as including this tag ensures that the documentation for these objects appears in the
appropriate part of the manual.  This page of the manual is linked to from all of the pages that exploit the reference
module's functionality to provide multiple methods for calculating the distance between two trajectory frames.  Any measure
that you implement should thus have one of these wrapper actions associated with it and the documentation for the measure
should be included in the wrapper code's source code file.

