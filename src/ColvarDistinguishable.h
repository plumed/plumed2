#ifndef __PLUMED_ColvarDistinguishable_h
#define __PLUMED_ColvarDistinguishable_h

#include <vector>
#include "Colvar.h"

namespace PLMD {

class ColvarDistinguishable : public Colvar {
private:
/// The indexes of the atoms in each colvar
  std::vector< std::vector<unsigned> > function_indexes;
/// The derivatives with respect to the atoms for a single thing in the list    
  std::vector<Vector> derivatives;
/// The virial with respect to the atoms for a single thing in the list
  Tensor virial;
protected:
/// Add some indexes to create another colvar to calculate
  void addIndexes( const unsigned& astart, const std::vector<unsigned>& new_indexes );
public:
  ColvarDistinguishable(const ActionOptions&);
/// Return the number of colvars
  unsigned getNumberOfColvars() const;
/// Routines used to transfer the derivatives for a single colvar onto the list of derivatives
  void mergeFunctions( const unsigned& vf, const unsigned& nf, const double& f, const double& df );
  void mergeFunctions( const std::string& valname, const unsigned& nf, const double& f, const double& df );
/// Interpret the use of the groups keyword
  void interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups );
/// Interpret the appearance of the atoms keyword
  void interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist );
/// Update the neighbour list
  void updateDynamicContent( const double& cutoff, std::vector<bool>& skips );
/// Calculate a single colvar
  double calcFunction( const unsigned& i );
/// Compute the value of the colvar
  virtual double compute( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial )=0;
};

inline
unsigned ColvarDistinguishable::getNumberOfColvars() const {
  return function_indexes.size();
}

inline
double ColvarDistinguishable::calcFunction( const unsigned& i ){
  assert( i<function_indexes.size() ); virial.clear();
  for(unsigned i=0;i<derivatives.size();++i) derivatives[i].clear();
  return compute( function_indexes[i], derivatives, virial );
}

inline
void ColvarDistinguishable::mergeFunctions( const unsigned& vf, const unsigned& nf, const double& f, const double& df ){
  const unsigned nat=getNumberOfAtoms();
  for(unsigned i=0;i<function_indexes[nf].size();++i){
      addAtomicDerivative( vf, function_indexes[nf][i], f, df*derivatives[i] );
  }
  addVirial( vf, f, df*virial );
}

inline
void ColvarDistinguishable::mergeFunctions( const std::string& valname, const unsigned& nf, const double& f, const double& df ){
  unsigned vf=getValueNumberForLabel( valname );
  const unsigned nat=getNumberOfAtoms();
  for(unsigned i=0;i<function_indexes[nf].size();++i){
      addAtomicDerivative( vf, function_indexes[nf][i], f, df*derivatives[i] );
  }
  addVirial( vf, f, df*virial );
}

}
#endif
