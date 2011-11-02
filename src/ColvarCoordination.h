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
void ColvarCoordination::mergeFunctions( const unsigned& vf, const unsigned& nf, const double& f, const double& df ){
  const unsigned nat=getNumberOfAtoms();
  
  addAtomicDerivative( vf, central[nf], f, df*derivatives[0] );
  unsigned nd=1;
  for(unsigned i=0;i<sphere[nf].skipto.size();i=sphere[nf].skipto[i]){
      addAtomicDerivative( vf, sphere[nf].index[i], f, df*derivatives[nd] );
      nd++;
  }
  addVirial( vf, f, df*virial );
}

inline
void ColvarCoordination::mergeFunctions( const std::string& valname, const unsigned& nf, const double& f, const double& df ){
  unsigned vf=getValueNumberForLabel( valname );
  const unsigned nat=getNumberOfAtoms();

  addAtomicDerivative( vf, central[nf], f, df*derivatives[0] );
  unsigned nd=1;
  for(unsigned i=0;i<sphere[nf].skipto.size();i=sphere[nf].skipto[i]){
      addAtomicDerivative( vf, sphere[nf].index[i], f, df*derivatives[nd] );
      nd++;
  }
  addVirial( vf, f, df*virial );
}

}
#endif

