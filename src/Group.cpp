#include "Group.h"

using namespace std;

namespace PLMD {

Group::Group(const ActionOptions&ao):
  ActionAtomistic(ao)
{
  allowKeyword("ATOMS"); forbidKeyword("STRIDE"); forbidKeyword("NL_CUTOFF");
}

Group::~Group(){
  plumed.getAtoms().removeGroup(getLabel());
}

void Group::readGroup(){
  int natoms=-1; unsigned ngrp=1;
  readActionAtomistic( natoms, ngrp ); std::vector<double> domain(2,0.0); 
  readActionWithExternalArguments( 3*getNumberOfAtoms()+9, domain );
  // The value is the number of atoms that satisfy the conditions defined by the group
  // be aware this could be a non-integer value due to the use of differentiable functions
  addValue("natoms", false, true );

  // Now set up the skipto list
  skipto.resize( getNumberOfAtoms() ); contributions.resize( natoms ); 
  derivatives.resize( natoms ); positions.resize( natoms );
  f.resize( natoms ); forces.resize( 3*natoms + 9 );
  for(unsigned i=0;i<getNumberOfAtoms();++i) skipto[i]=i+1;
}

void Group::interpretAtomsKeyword( const std::vector<std::vector<unsigned> >& flist ){
  if( flist.size()!=1 ) error("cannot create multiple groups in a single line");

  log.printf("  created from the following list of atoms : " );
  for(unsigned j=0;j<flist[0].size();++j) log.printf("%s ", plumed.getAtoms().interpretIndex( flist[0][j] ).c_str() );
  log.printf("\n");
}

void Group::updateDynamicContent( const double& cutoff, std::vector<bool>& skips ){
  assert( skips.size()==getNumberOfAtoms() );

  positions.resize( getNumberOfAtoms() );
  contributions.resize( getNumberOfAtoms() );
  derivatives.resize( getNumberOfAtoms() );
  for(unsigned i=0;i<getNumberOfAtoms();++i){
     positions[i]=getPositions(i); derivatives[i].clear();
  }  
  Tensor v;
  double nvla=compute( positions, contributions, derivatives, v );

  unsigned n=0, nactive=0;
  for(unsigned i=1;i<getNumberOfAtoms();++i){
      assert( contributions[i]<=1.0 );   // Should be a number between 1 and 0 always
      if( contributions[i]<cutoff ){
          skips[i]=true;
      } else if( !skips[i] ){
         skipto[n]=i; n=i; nactive++;
      }
  }   
  skipto[n]=getNumberOfAtoms();
  // Resize everything
  positions.resize( nactive ); contributions.resize( nactive ); derivatives.resize( nactive );
}

void Group::calculate(){
  calculateAtomisticActions();

  // Get the positions that are not being skipped due to neighbour lists
  unsigned n=0;
  for(unsigned i=0;i<skipto.size();i=skipto[i]){
     positions[n]=getPositions(i); derivatives[n].clear(); n++;
  }
  assert( n==positions.size() );

  // Now do all the computations
  Tensor v; n=0;
  double nvla=compute( positions, contributions, derivatives, v );
  for(unsigned i=0;i<skipto.size();i=skipto[i] ){
      addAtomicDerivative( 0, i, contributions[n], derivatives[n] );
      n++;
  }
  assert( n==contributions.size() ); addVirial( 0, nvla, v );
  setValue( 0, nvla, 1.0 );
}

void Group::apply(){
  const unsigned nat=f.size();
  for(unsigned i=0;i<nat;i++){
    f[i][0]=0.0; f[i][1]=0.0; f[i][2]=0.0;
  }

  Tensor v; v.clear();
  for(int i=0;i<getNumberOfValues();++i){
    if( getForces( i, forces ) ){
       for(unsigned j=0;j<nat;++j){
          f[j][0]+=forces[3*j+0];
          f[j][1]+=forces[3*j+1];
          f[j][2]+=forces[3*j+2];
       }
       v(0,0)+=forces[3*nat+0];
       v(0,1)+=forces[3*nat+1];
       v(0,2)+=forces[3*nat+2];
       v(1,0)+=forces[3*nat+3];
       v(1,1)+=forces[3*nat+4];
       v(1,2)+=forces[3*nat+5];
       v(2,0)+=forces[3*nat+6];
       v(2,1)+=forces[3*nat+7];
       v(2,2)+=forces[3*nat+8];
    }
  } 
  applyForces( f, v );
}

}
