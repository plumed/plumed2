#ifndef __PLUMED_ColvarCoordination_h
#define __PLUMED_ColvarCoordination_h

#include <string>
#include <cmath>
#include <cassert>
#include <vector>
#include "Colvar.h"

namespace PLMD {

class csphere {
friend class ColvarCoordination;
private:
  std::vector<unsigned> index;
  std::vector<unsigned> skipto;
public:
  void resize( const unsigned sz ){ index.resize(sz); skipto.resize(sz); }
};

class ColvarCoordination : public Colvar {
private:
/// The list of central atoms
  std::vector<unsigned> central;
/// The list of coordination spheres
  std::vector<csphere> sphere;
/// The list of atom positions
  std::vector<Vector> neighbours;
/// The derivatives with respect to the vectors connecting the central atom to the coordination sphere
  std::vector<Vector> derivatives;
/// The derivatives with respect to the cell box
  Tensor virial;
public:
  ColvarCoordination(const ActionOptions&);
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
/// Compute a single colvar
  virtual double compute( const unsigned& nvecs, const Vector& cpos, const std::vector<Vector>& npos, std::vector<Vector>& derivatives, Tensor& vir )=0;
};

inline
unsigned ColvarCoordination::getNumberOfColvars() const {
  return central.size();
}

inline
double ColvarCoordination::calcFunction( const unsigned& i ){
  assert( i<central.size() ); unsigned k=0;
  derivatives[k][0]=derivatives[k][1]=derivatives[k][2]=0.0;
  for(unsigned j=0;j<sphere[i].skipto.size();j=sphere[i].skipto[j]){
     neighbours[k]= getPositions( sphere[i].index[j] ); k++;
     derivatives[k][0]=derivatives[k][1]=derivatives[k][2]=0.0;
  }
  virial.clear(); 
  return compute( k, getPositions( central[i] ), neighbours, derivatives, virial );
}

inline
void ColvarCoordination::mergeFunctions( const unsigned& vf, const unsigned& nf, const double& df ){
  const unsigned nat=getNumberOfAtoms();
  
  addDerivative( vf, 3*central[nf] + 0 , df*derivatives[0][0] );
  addDerivative( vf, 3*central[nf] + 1 , df*derivatives[0][1] );
  addDerivative( vf, 3*central[nf] + 2 , df*derivatives[0][2] );
  unsigned nd=1;
  for(unsigned i=0;i<sphere[nf].skipto.size();i=sphere[nf].skipto[i]){
      addDerivative( vf, 3*sphere[nf].index[i] + 0 ,  df*derivatives[nd][0] );
      addDerivative( vf, 3*sphere[nf].index[i] + 1 ,  df*derivatives[nd][1] );
      addDerivative( vf, 3*sphere[nf].index[i] + 2 ,  df*derivatives[nd][2] );
      nd++;
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
void ColvarCoordination::mergeFunctions( const std::string& valname, const unsigned& nf, const double& df ){
  unsigned vf=getValueNumberForLabel( valname );
  const unsigned nat=getNumberOfAtoms();

  addDerivative( vf, 3*central[nf] + 0 , df*derivatives[0][0] );
  addDerivative( vf, 3*central[nf] + 1 , df*derivatives[0][1] );
  addDerivative( vf, 3*central[nf] + 2 , df*derivatives[0][2] );
  unsigned nd=1;
  for(unsigned i=0;i<sphere[nf].skipto.size();i=sphere[nf].skipto[i]){
      addDerivative( vf, 3*sphere[nf].index[i] + 0 ,  df*derivatives[nd][0] );
      addDerivative( vf, 3*sphere[nf].index[i] + 1 ,  df*derivatives[nd][1] );
      addDerivative( vf, 3*sphere[nf].index[i] + 2 ,  df*derivatives[nd][2] );
      nd++;
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

