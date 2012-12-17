/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "DRMSD.h"
#include "PDB.h"
#include "Pbc.h"
#include <cmath>

using namespace std;
namespace PLMD{

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

  plumed_assert(positions.size()==natoms && derivatives.size()==natoms );

  Vector distance;
  double drmsd=0.; virial.zero();
  for(unsigned i=0;i<derivatives.size();++i) derivatives[i].zero();
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
    virial -=  ( diff / len ) * Tensor(distance,distance); 
  }

  double npairs = static_cast<double>(targets.size());

  drmsd=sqrt( drmsd / npairs );
    
  double idrmsd=1.0/( drmsd * npairs );

  virial *= idrmsd;

  for(unsigned i=0;i<natoms;++i){derivatives[i] *= idrmsd;}

  return drmsd;
}


}
