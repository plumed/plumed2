#ifndef __PLUMED_MultiColvar_h
#define __PLUMED_MultiColvar_h

#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include "ActionWithDistribution.h"
#include <vector>

#define PLUMED_MULTICOLVAR_INIT(ao) Action(ao),MultiColvar(ao)

namespace PLMD {

class MultiColvar;

class AtomList {
private:
  MultiColvar* mcolvar;
  std::vector<unsigned> all_atoms;
  std::vector<unsigned> current_atoms;
public:
  AtomList( MultiColvar* c );
  void addAtom( const unsigned n ); 
  void clear();
  void getAtoms(std::vector<Vector>& positions ) const ;
  unsigned getNumberOfAtoms() const;
  unsigned getAtomNumber( const unsigned& i ) const ;
};

inline
unsigned AtomList::getNumberOfAtoms() const {
  return current_atoms.size();
}

inline
unsigned AtomList::getAtomNumber( const unsigned& i ) const {
  plumed_massert(i<current_atoms.size(),"There are not enough atoms");
  return all_atoms[ current_atoms[i] ];
}

/// Action representing a collective variable
class MultiColvar :
  public ActionAtomistic,
  public ActionWithValue,
  public ActionWithDistribution
  {
friend class AtomList;
private:
  bool usepbc;
  bool readatoms;
  std::vector<unsigned> index_translator;
  std::vector<AtomNumber> real_atoms;
  std::vector<AtomList> ColvarAtoms;
protected:
/// Read the keyword ATOMS
  void readAtomsKeyword( int& natoms );
/// Update the atoms request
  void requestAtoms();
/// Get the separation between a pair of vectors
  Vector getSeparation( const Vector& vec1, const Vector& vec2 ) const ;
/// Get the position of a particular atom in the array
  Vector getAtomPosition( const unsigned j ) const ;
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
/// Merge the derivatives 
  void mergeDerivatives( const unsigned j, Value* value_in, Value* value_out );
/// Calcualte the colvar
  void calculateThisFunction( const unsigned& j, Value* value_in );
/// And a virtual function which actually computes the colvar
  virtual double compute( const std::vector<Vector>& pos, std::vector<Vector>& deriv, Tensor& virial )=0;  
};

inline
unsigned MultiColvar::getNumberOfFunctionsInDistribution(){
  return ColvarAtoms.size();
}

inline
unsigned MultiColvar::getThisFunctionsNumberOfDerivatives( const unsigned& j ){
  return 3*ColvarAtoms[j].getNumberOfAtoms() + 9;
}

inline
Vector MultiColvar::getAtomPosition( const unsigned j ) const {
  plumed_massert(j<real_atoms.size(),"Atoms has gone out of bounds");
  unsigned atom=index_translator[j];
  plumed_massert(atom>=0,"Atom has currently not been requested");
  return getPosition(atom);
}

}

#endif
