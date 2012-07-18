#include "DRMSD.h"
#include "PDB.h"
#include <cassert>
#include <cmath>

using namespace std;
using namespace PLMD;

void DRMSD::setFromPDB( const double& bondc, const PDB&pdb ){
  setReference(bondc, pdb.getPositions());
}

void DRMSD::clear(){
  targets.resize(0,0);
}

void DRMSD::setReference( const double& bondc, const vector<Vector> & reference ){
  natoms=reference.size(); npairs=0;
  targets.resize( natoms, natoms );
  for(unsigned i=1;i<natoms;++i){
     for(unsigned j=0;j<i;++j){
         targets(i,j)=delta(reference[i],reference[j]).modulo();
         if( targets(i,j)<bondc ){
             targets(i,j)=-4.0;
         } else {
             npairs++;
         }
     }
  }
}

double DRMSD::calculate(const std::vector<Vector> & positions,std::vector<Vector> &derivatives, Tensor& virial) const {
  assert(positions.size()==natoms && derivatives.size()==natoms );

  Vector distance; double len,diff,drmsd=0; 
  for(unsigned i=1;i<natoms;++i){
     for(unsigned j=0;j<i;++j){
        if( targets(i,j)>0 ){
            distance=delta( positions[i] , positions[j] ); 
            len=distance.modulo();  
            diff=len-targets(i,j); drmsd+=diff*diff; 
            derivatives[i]+=-( diff / len ) * distance;
            derivatives[j]+= ( diff / len ) * distance;
            virial=virial-( diff / len ) * Tensor(distance,distance);
        }
     }
  }
  drmsd=sqrt( drmsd / npairs ); 
     
  double idrmsd=1.0/( drmsd * npairs ); virial*=idrmsd;
  for(unsigned i=0;i<natoms;++i){
      derivatives[i]=idrmsd*derivatives[i]; 
  }

  return drmsd;
}

double DRMSD::calculate(const std::vector<Vector> & positions, const Pbc& pbc, std::vector<Vector> &derivatives, Tensor& virial) const {
  assert(positions.size()==natoms && derivatives.size()==natoms );

  Vector distance; double len,diff,drmsd=0;
  for(unsigned i=1;i<natoms;++i){
     for(unsigned j=0;j<i;++j){
        if( targets(i,j)>0 ){
            distance=pbc.distance( positions[i] , positions[j] ); 
            len=distance.modulo();  
            diff=len-targets(i,j); drmsd+=diff*diff;
            derivatives[i]+=-( diff / len ) * distance;
            derivatives[j]+= ( diff / len ) * distance;
            virial=virial-( diff / len ) * Tensor(distance,distance); 
        }
     }
  }
  drmsd=sqrt( drmsd / npairs );
    
  double idrmsd=1.0/( drmsd * npairs ); virial*=idrmsd;
  for(unsigned i=0;i<natoms;++i){
      derivatives[i]=idrmsd*derivatives[i]; 
  }

  return drmsd;
}


