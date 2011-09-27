#ifndef __PLUMED_Colvar_h
#define __PLUMED_Colvar_h

#include "ActionAtomistic.h"
#include "ActionWithValue.h"
#include "SwitchingFunction.h"
#include "HistogramBead.h"
#include <vector>

namespace PLMD {

/// Action representing a collective variable
class Colvar : public ActionAtomistic {
private:
/// The indexes of the atoms in each colvar
  std::vector< std::vector<unsigned> > function_indexes;
/// Makes the neighbour list work
  std::vector<unsigned> skipto;
/// The derivatives with respect to the atoms for a single thing in the list    
  std::vector<Vector> derivatives;
/// The virial with respect to the atoms for a single thing in the list
  Tensor virial;
/// The forces on the atoms and on the virial
  std::vector<double> forces;
/// The forces on the atoms
  std::vector<Vector> f;
/// What kind of calculation are we doing
  bool doall, domin, domax, dototal, domean, dolt, domt, dohist;
/// The beta parameter for caclulating the minimum
  double beta;
/// Switching function for less than and more than    
  SwitchingFunction ltswitch, mtswitch;
/// The beads for the histograms and within
  std::vector<HistogramBead> histogram;
/// This is the vector we store the histogram beads inside
  std::vector<double> hvalues;
/// Routines used to transfer the derivatives for a single colvar onto the list of derivatives
  void mergeFunctions( const unsigned& vf, const unsigned& nf, const double& df );
  void mergeFunctions( const std::string& valname, const unsigned& nf, const double& df );
protected:
  void readActionColvar( int natoms, const std::vector<double>& domain );
public:
  Colvar(const ActionOptions&);
  ~Colvar(){};
  void interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups );
  void interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist );
  void updateNeighbourList( const double& cutoff, std::vector<bool>& skips );
  void calculate();
  void apply();
  virtual double calcFunction( const std::vector<unsigned>& indexes, std::vector<Vector>& derivatives, Tensor& virial )=0; 
};

inline
void Colvar::mergeFunctions( const unsigned& vf, const unsigned& nf, const double& df ){
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
void Colvar::mergeFunctions( const std::string& valname, const unsigned& nf, const double& df ){
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
