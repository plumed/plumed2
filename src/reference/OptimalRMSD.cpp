/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013,2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

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
#include "MetricRegister.h"
#include "RMSDBase.h"
#include "tools/Matrix.h"
#include "tools/RMSD.h"

namespace PLMD{

class OptimalRMSD : public RMSDBase {
private:
  bool fast;
  RMSD myrmsd;
public:
  OptimalRMSD(const ReferenceConfigurationOptions& ro);
  void read( const PDB& );
  double calc( const std::vector<Vector>& pos, const bool& squared );
};

PLUMED_REGISTER_METRIC(OptimalRMSD,"OPTIMAL")

OptimalRMSD::OptimalRMSD(const ReferenceConfigurationOptions& ro ):
ReferenceConfiguration(ro),
RMSDBase(ro)
{
  fast=ro.usingFastOption();
}

void OptimalRMSD::read( const PDB& pdb ){
  readReference( pdb ); 
}

double OptimalRMSD::calc( const std::vector<Vector>& pos, const bool& squared ){
  if( fast ){
     if( getAlign()==getDisplace() ) return myrmsd.optimalAlignment<false,true>(getAlign(),getDisplace(),pos,getReferencePositions(),atom_ders,squared); 
     return myrmsd.optimalAlignment<false,false>(getAlign(),getDisplace(),pos,getReferencePositions(),atom_ders,squared);
  } else {
     if( getAlign()==getDisplace() ) return myrmsd.optimalAlignment<true,true>(getAlign(),getDisplace(),pos,getReferencePositions(),atom_ders,squared);
     return myrmsd.optimalAlignment<true,false>(getAlign(),getDisplace(),pos,getReferencePositions(),atom_ders,squared);
  }
}

}
