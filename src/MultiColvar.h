#ifndef __PLUMED_MultiColvar_h
#define __PLUMED_MultiColvar_h

#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include "ActionWithDistribution.h"
#include "NeighbourList.h"
#include <vector>

#define PLUMED_MULTICOLVAR_INIT(ao) Action(ao),MultiColvar(ao)

namespace PLMD {

class MultiColvar;

/// Action representing a collective variable
class MultiColvar :
  public ActionAtomistic,
  public ActionWithValue,
  public ActionWithDistribution
  {
private:
  bool usepbc;
  bool readatoms;
/// The list of all the atoms involved in the colvar
  DynamicList all_atoms;
/// The lists of the atoms involved in each of the individual colvars
/// note these refer to the atoms in all_atoms
  std::vector<DynamicList> colvar_atoms;
/// The stuff to control the frequency of the neighbour list update
  bool setupList;
  int updateFreq;
  unsigned lastUpdate;
  bool reduceAtNextStep;
/// The neighbour list we use the same one for every colvar
  NeighbourList<Vector,Pbc> nlist;
/// Read in ATOMS keyword
  void readAtomsKeyword( int& natoms );
/// Read in the various GROUP keywords
  void readGroupsKeyword( int& natoms );
/// Read in the various SPECIES keywords
  void readSpeciesKeyword( int& natoms );
/// Retrieve the positions of the atoms
  void retrieveAtoms( const unsigned& j, std::vector<Vector>& pos );
/// Update the atoms request
  void requestAtoms();
protected:
/// Read in all the keywords that can be used to define atoms
  void readAtoms( int& natoms );
/// Create a neighbour list that is composed of these atoms
  void createNeighbourList( std::vector<std::pair<unsigned,unsigned> >& pairs );
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
public:
  MultiColvar(const ActionOptions&);
  ~MultiColvar(){};
  static void registerKeywords( Keywords& keys );
  static void useNeighbourList( const std::string& style, Keywords& keys );
/// If it is a neighbour list update step make sure we are requesting all atoms
  void prepare();
/// Apply the forces on the values
  void apply();
/// Return the number of Colvars this is calculating
  unsigned getNumberOfFunctionsInDistribution();  
/// Return the number of derivatives for a given colvar
  unsigned getThisFunctionsNumberOfDerivatives( const unsigned& j ); 
/// Merge the derivatives 
  void mergeDerivatives( const unsigned j, Value* value_in, Value* value_out );
/// Calcualte the colvar
  void calculateThisFunction( const unsigned& j, Value* value_in );
/// And a virtual function which actually computes the colvar
  virtual double compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial )=0;  
};

inline
unsigned MultiColvar::getNumberOfFunctionsInDistribution(){
  return colvar_atoms.size();
}

inline
unsigned MultiColvar::getThisFunctionsNumberOfDerivatives( const unsigned& j ){
  return 3*colvar_atoms[j].getNumberActive() + 9;
}

}

#endif
