#include "DRMSD.h"
#include "PDB.h"
#include <cassert>
#include <cmath>

using namespace std;
using namespace PLMD;

void DRMSD::setFromPDB(const PDB&pdb, double lbound, double ubound){
  setReference(pdb.getPositions(), lbound, ubound);
}

void DRMSD::clear(){
  targets.clear();
}

void DRMSD::setReference(const vector<Vector> & reference, double lbound, double ubound)
{
 natoms = reference.size();
 for(unsigned i=0;i<natoms-1;++i){
  for(unsigned j=i+1;j<natoms;++j){
   double distance = delta(reference[i],reference[j]).modulo();
   if(distance > lbound && distance < ubound){ 
     targets[make_pair(i,j)] = distance; 
   }
  }
 }
}

double DRMSD::calculate(const std::vector<Vector> & positions,
                        std::vector<Vector> &derivatives, Tensor& virial) const {

 Pbc fake_pbc;
 return DRMSD::calculate(positions,fake_pbc,derivatives,virial,false);
}

double DRMSD::calculate(const std::vector<Vector> & positions, const Pbc& pbc,
                        std::vector<Vector> &derivatives, Tensor& virial, bool do_pbc) const {

  assert(positions.size()==natoms && derivatives.size()==natoms );

  Vector distance;
  double drmsd=0.;
  for(map< pair <unsigned,unsigned> , double>::const_iterator it=targets.begin();it!=targets.end();++it){
    unsigned i=it->first.first;
    unsigned j=it->first.second;
    if(do_pbc){distance=pbc.distance( positions[i] , positions[j] );}
    else{distance=delta( positions[i] , positions[j] );}
    double len = distance.modulo();  
    double diff = len - it->second;
    drmsd += diff * diff;
    derivatives[i] += -( diff / len ) * distance;
    derivatives[j] +=  ( diff / len ) * distance;
    virial -= ( diff / len ) * Tensor(distance,distance); 
  }

  double npairs = static_cast<double>(targets.size());

  drmsd=sqrt( drmsd / npairs );
    
  double idrmsd=1.0/( drmsd * npairs );
  virial*=idrmsd;
  for(unsigned i=0;i<natoms;++i){derivatives[i] *= idrmsd;}

  return drmsd;
}


