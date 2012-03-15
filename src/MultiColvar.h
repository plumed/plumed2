#ifndef __PLUMED_MultiColvar_h
#define __PLUMED_MultiColvar_h

#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include "ActionWithDistribution.h"
#include <vector>

#define PLUMED_MULTICOLVAR_INIT(ao) Action(ao),MultiColvar(ao)

namespace PLMD {

class MultiColvar;

/// Action for multiple collective coordinates and distributions
/// of collective variables.

/**
As you are no doubt aware within plumed 2 you can calculate multiple 
instances of a collective coorinate from a single line in the input file.
One can then calculate functions such as the minimum, number less than,...
from the resulting distribution of collective variables.  To create these kinds
of collective variables we use the functionality implmented in PLMD::MultiColvar.
In fact in writing a single PLMD::MultiColvar you actually write many CVs at once
as the minimum of a distribution, the number less than and so on come with no 
additional effort on your part.

To better understand how to go about writing a new MultiColvar examine the interior
of one of the existing MultiColvars, e.g. MultiColvarDistance.cpp or MultiColvarCoordination.cpp.
In fact a good way to start is to copy one of these files as certain features (e.g. the fact that you
have to include MultiColvar.h and ActionRegister.h and that all your code should be inside the namespace
PLMD) never change. 

\section docs Creating documentation  

The first thing you will have to change is the documentation.  As discussed on the \ref usingDoxygen page
of this manual the documentation is created using Doxygen.  You are implementing a cv so your PLMEDOC line
should read:

\verbatim
//+PLUMEDOC COLVAR MYCVKEYWORD 
\endverbatim 

Your documentation should contain a description of what your CV calculates and some examples.  You do not
need to write a description of the input syntax as that will be generated automatically.  For more information
on how to write documentation go to \ref usingDoxygen.

\section class Creating the class

The first step in writing the executable code for your CV is to create a class that calculates your CV.  The
declaration for your class should appear after the documentation and will look something like:

\verbatim
class MyNewMultiColvar : public MultiColvar {
private:
   // Declare all the variables you need here  
public:
  static void registerKeywords( Keywords& keys );
  MyNewMultiColvar(const ActionOptions&);
  virtual double compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial );
};

PLUMED_REGISTER_ACTION(MyNewMultiColvar,"MYCVKEYWORD")
\endverbatim

This new class (MyNewMultiColvar) inherits from MultiColvar and so contains much of the functionality we 
require already.  Furthermore, by calling PLUMED_REGISTER_ACTION we have ensured that whenever the keyword
MYCVKEYWORD is found in the input file an object of type MyNewMultiColvar will be generated and hence that your 
new CV will be calculated wherever it is required.

The three functions that are defined in the above class (registerKeywords, the constructor and compute) are mandatory.  
Without these functions the code will not compile and run.  Writing your new CV is simply a matter of writing these
three subroutines.

\section register RegisterKeywords

RegisterKeywords is the routine that is used by plumed to create the remainder of the documentation.  As much of this
documentation is created inside the MultiColvar class itself the first line of your new registerKeywords routine must
read:

\verbatim
MultiColvar::registerKeywords( keys )
\endverbatim

as this creates all the documentation for the features that are part of all PLMD::MultiColvar objects.  To see how to create 
documentation that describes the keywords that are specific to your particular CV please read \ref usingDoxygen.

\par Reading the atoms

Creating the lists of atoms involved in each of the colvars in a PLMD::MultiColvar is quite involved.  What is more this is
an area where we feel it is important to maintain some consistency in the input.  For these reasons one must decide at keyword
registration time what the most appropriate way of reading the atoms is for your new CV from the following list of options:

<table align=center frame=void width=95%% cellpadding=5%%>
<tr>
<td width=5%> ATOMS </td> <td> The atoms keyword specifies that one collective coordinate is to be calculated for each set of atoms specified.  
                               Hence, for MultiColvarDistance the command MULTIDISTANCE ATOMS1=1,2 ATOMS2=2,3 ATOMS3=3,4 specifies 
                               that three distances should be calculated. </td>
</tr> <tr>
<td width=5%> GROUP GROUPA GROUPB </td> <td> The GROUP keyword is used for quantities such as distances and angles.  A single GROUP specifies that
                                             a CV should be calculated for each distinct set of atoms that can be made from the group.  Hence,
                                             MUTIDISTANCE GROUP=1,2,3 specifies that three distances should be calculated (1,2), (1,3) and (2,3).
                                             If there is a GROUPA and GROUPB one CV is calculated for each set of atoms that includes at least one
                                             member from each group.  Thus MULTIDISTANCE GROUPA=1 GROUPB=2,3 calculates two distance (1,2) and (1,3). </td>
</tr> <tr>
<td width=5%> SPECIES SPECIESA SPECIESB </td> <td> The SPECIES keywords is used for quantities like coordination numbers.  The way this works is
                                                   best explained using an example.  Imagine a user working on NaCl wishes to calculate the average
                                                   coordination number of the sodium ions in his/her system.  To do this he/she uses COORDINATIONNUMBER SPECIES=1-100
                                                   which tells plumed that atoms 1-100 are the sodium ions.  To calculate the average coordination number
                                                   plumed calculates 100 coordination numbers (one for each sodium) and averages.  Obviously, each of these coordination
                                                   numbers involves the full set of 100 sodium atoms.  By contrast if the user wanted to calculate the coordination number
                                                   of Na with Cl he/she would do COORDINATIONNUMBER SPECIESA=1-100 SPECIESB=101-200, where obviously 101-200 are the 
                                                   chlorine ions.  Each of these 100 heteronuclear coordination numbers involves the full set of atoms speciefied using
                                                   SPECIESB and one of the ions speciefied by SPECIESA. </td>
</tr>
</table>

You may use one or more of these particular options by doing:

\verbatim
keys.use("ATOMS");
\endverbatim

for all the keywords in the first cell of the particular record that you want to use from the table above.

\par Parallelism

The most optimal way to parallelize your CV will depend on what it is your CV calculates.  In fact it will probably even
depend on what the user puts in the input file.  If you know what you are doing feel free to parallelize your CV 
however you feel is best. However, the simplest way to parallize a long list of CVs (a multi-colvar) and functions thereof
is to parallelize over CVs.  Even if you are unsure of what you are doing with MPI you can turn this form of parallelism on
by calling PLMD::ActionWithDistribution::autoParallelize from within your registerKeywords functions.     

\section label4 The constructor

The constructor reads the plumed input files and as such it must:

- Read all the keywords that are specific to your new CV 
- Read in the atoms involved in the CV
- And finally check that readin was succesful

The first of these tasks is dealt with in detail on the following page: \usingDoxygen.  Furthermore, to do the second and third 
tasks on this list one simply calls PLMD::MultiColvar::readAtoms and PLMD::Action::checkRead.  (you must have decided which of the 
keywords described above you are using to read the atoms however for readAtoms to work).  

\section label5 Compute

Compute is the routine that actually calculates the colvar.  It takes as input an array of positions and returns the derivatives with
repect to the atomic positions, the virial and the value of the colvar (this is the double that the routine returns).  When calculating
multiple colvars using a MultiColvar this routine will be called more than once with different sets of atomic positions.  The input std::vector
of positions has the atomic positions ordered in the way they were specified to the ATOMS keyword.  If SPECIES is used the central atom is 
in the first position of the position vector and the atoms in the coordination sphere are in the remaining elements.  

*/

class MultiColvar :
  public ActionAtomistic,
  public ActionWithValue,
  public ActionWithDistribution
  {
private:
  bool usepbc;
  bool readatoms;
  bool needsCentralAtomPosition;
/// The list of all the atoms involved in the colvar
  DynamicList<AtomNumber> all_atoms;
/// The lists of the atoms involved in each of the individual colvars
/// note these refer to the atoms in all_atoms
  std::vector< DynamicList<unsigned> > colvar_atoms;
/// Used to stop calculating during neighbor list update
  bool stopcondition;
/// Used to make sure we update the correct atoms during neighbor list update
  unsigned current;
/// Read in ATOMS keyword
  void readAtomsKeyword( int& natoms );
/// Read in the various GROUP keywords
  void readGroupsKeyword( int& natoms );
/// Read in the various SPECIES keywords
  void readSpeciesKeyword( int& natoms );
/// Update the atoms request
  void requestAtoms();
protected:
/// Read in all the keywords that can be used to define atoms
  void readAtoms( int& natoms );
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Stop a calculation during neighbor list update steps
  void stopCalculatingThisCV();
/// Update the list of atoms after the neighbor list step
  void removeAtomRequest( const unsigned& aa );
/// Do we use pbc to calculate this quantity
  bool usesPbc() const ;
public:
  MultiColvar(const ActionOptions&);
  ~MultiColvar(){};
  static void registerKeywords( Keywords& keys );
/// Apply the forces on the values
  void apply();
/// Return the number of Colvars this is calculating
  unsigned getNumberOfFunctionsInDistribution();  
/// Return the number of derivatives for a given colvar
  unsigned getThisFunctionsNumberOfDerivatives( const unsigned& j ); 
/// Expand the lists of atoms so we get them all when it is time to update the neighbor lists
  void prepareForNeighborListUpdate();
/// Contract the lists of atoms once we have finished updating the neighbor lists so we only get
/// the required subset of atoms 
  void completeNeighborListUpdate();
/// Merge the derivatives 
  void mergeDerivatives( const unsigned j, Value* value_in, Value* value_out );
/// Turn of atom requests when this colvar is deactivated cos its small
  void deactivate( const unsigned j );
/// Calcualte the colvar
  void calculateThisFunction( const unsigned& j, Value* value_in, std::vector<Value>& aux );
/// And a virtual function which actually computes the colvar
  virtual double compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial )=0; 
/// A virtual routine to get the position of the central atom - used for things like cv gradient
  virtual void getCentralAtom( const std::vector<Vector>& pos, std::vector<Value>& pos); 
};

inline
unsigned MultiColvar::getNumberOfFunctionsInDistribution(){
  return colvar_atoms.size();
}

inline
void MultiColvar::deactivate( const unsigned j ){
  colvar_atoms[j].deactivateAll();
}

inline
unsigned MultiColvar::getThisFunctionsNumberOfDerivatives( const unsigned& j ){
  return 3*colvar_atoms[j].getNumberActive() + 9;
}

inline
void MultiColvar::stopCalculatingThisCV(){
  plumed_massert(isTimeForNeighborListUpdate(),"found stop but not during neighbor list step");
  stopcondition=true;
}

inline
void MultiColvar::removeAtomRequest( const unsigned& i ){
  plumed_massert(isTimeForNeighborListUpdate(),"found removeAtomRequest but not during neighbor list step");
  colvar_atoms[current].deactivate(i); 
}

inline
bool MultiColvar::usesPbc() const {
  return usepbc;
}

}

#endif
