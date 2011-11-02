#ifndef __PLUMED_Group_h
#define __PLUMED_Group_h

#include "ActionAtomistic.h"
#include "PlumedMain.h"

namespace PLMD{

class Group : public ActionAtomistic {
  std::vector<double> forces;
  std::vector<unsigned> skipto;
  std::vector<double> contributions;
  std::vector<Vector> f, positions, derivatives;
protected:
/// Read stuff related to group 
  void readGroup();
public:
  Group(const ActionOptions&ao);
  ~Group();
  void interpretGroupsKeyword( const unsigned& natoms, const std::string& atomGroupName, const std::vector<std::vector<unsigned> >& groups ){ assert(false); }
  void interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist );

/// Updates the atom selection in the group
  void updateDynamicContent( const double& cutoff, std::vector<bool>& skips );
/// Get the values of all the functions and the derivatives
  double getGroupData( std::vector<double>& group_f, std::vector<Vector>& group_df, Tensor& group_vir ) const ;
/// Updates the group of atoms involved
  virtual double compute( const std::vector<Vector>& positions, std::vector<double>& contributions, std::vector<Vector>& derivatives, Tensor& virial )=0;
/// Calcualte something
  void calculate();
/// Apply the derivatives
  void apply();
};

inline
double Group::getGroupData( std::vector<double>& group_f, std::vector<Vector>& group_df, Tensor& group_vir ) const {
  assert( group_f.size()==getNumberOfAtoms() && group_df.size()==getNumberOfAtoms() );

  // Set everything to zero
  for(unsigned i=0;i<getNumberOfAtoms();++i){ group_f[i]=0; group_df[i].clear(); }

  // Now get the ones we are not skipping
  unsigned n=0;
  for(unsigned i=0;i<getNumberOfAtoms();i=skipto[i]){
     group_f[i]=contributions[n]; group_df[i]=derivatives[n]; n++;
  }
  assert( n==derivatives.size() );

  unsigned nat=getNumberOfAtoms();
  group_vir(0,0)=getDerivative( 0, 3*nat + 0 );
  group_vir(0,1)=getDerivative( 0, 3*nat + 1 );
  group_vir(0,2)=getDerivative( 0, 3*nat + 2 );
  group_vir(1,0)=getDerivative( 0, 3*nat + 3 );
  group_vir(1,1)=getDerivative( 0, 3*nat + 4 );
  group_vir(1,2)=getDerivative( 0, 3*nat + 5 );
  group_vir(2,0)=getDerivative( 0, 3*nat + 6 );
  group_vir(2,1)=getDerivative( 0, 3*nat + 7 );
  group_vir(2,2)=getDerivative( 0, 3*nat + 8 );
  return getValue(0);
}

}

#endif
