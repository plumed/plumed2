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
  void mergeFunctions( const unsigned& vf, const unsigned& nf, const double& df );
  void mergeFunctions( const std::string& valname, const unsigned& nf, const double& df );
/// Interpret the use of the groups keyword
  void interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups );
/// Interpret the appearance of the atoms keyword
  void interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist );
/// Update the neighbour list
  void updateNeighbourList( const double& cutoff, std::vector<bool>& skips );
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
void ColvarDistinguishable::mergeFunctions( const unsigned& vf, const unsigned& nf, const double& df ){
  const unsigned nat=getNumberOfAtoms();
  for(unsigned i=0;i<function_indexes[nf].size();++i){
      addDerivative( vf, 3*function_indexes[nf][i] + 0 , df*derivatives[i][0] );
      addDerivative( vf, 3*function_indexes[nf][i] + 1 , df*derivatives[i][1] );
      addDerivative( vf, 3*function_indexes[nf][i] + 2 , df*derivatives[i][2] );
  }
  addDerivative( vf, 3*nat + 0, df*virial(0,0) );
  addDerivative( vf, 3*nat + 1, df*virial(0,1) );
  addDerivative( vf, 3*nat + 2, df*virial(0,2) );
  addDerivative( vf, 3*nat + 3, df*virial(1,0) );
  addDerivative( vf, 3*nat + 4, df*virial(1,1) );
  addDerivative( vf, 3*nat + 5, df*virial(1,2) );
  addDerivative( vf, 3*nat + 6, df*virial(2,0) );
  addDerivative( vf, 3*nat + 7, df*virial(2,1) );
  addDerivative( vf, 3*nat + 8, df*virial(2,2) );
}

inline
void ColvarDistinguishable::mergeFunctions( const std::string& valname, const unsigned& nf, const double& df ){
  unsigned vf=getValueNumberForLabel( valname );
  const unsigned nat=getNumberOfAtoms();
  for(unsigned i=0;i<function_indexes[nf].size();++i){
      addDerivative( vf, 3*function_indexes[nf][i] + 0 , df*derivatives[i][0] );
      addDerivative( vf, 3*function_indexes[nf][i] + 1 , df*derivatives[i][1] );
      addDerivative( vf, 3*function_indexes[nf][i] + 2 , df*derivatives[i][2] );
  }
  addDerivative( vf, 3*nat + 0, df*virial(0,0) );
  addDerivative( vf, 3*nat + 1, df*virial(0,1) );
  addDerivative( vf, 3*nat + 2, df*virial(0,2) );
  addDerivative( vf, 3*nat + 3, df*virial(1,0) );
  addDerivative( vf, 3*nat + 4, df*virial(1,1) );
  addDerivative( vf, 3*nat + 5, df*virial(1,2) );
  addDerivative( vf, 3*nat + 6, df*virial(2,0) );
  addDerivative( vf, 3*nat + 7, df*virial(2,1) );
  addDerivative( vf, 3*nat + 8, df*virial(2,2) );
}

}
#endif
