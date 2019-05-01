/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

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
#include "ActionVolume.h"

namespace PLMD {
namespace multicolvar {

void ActionVolume::registerKeywords( Keywords& keys ) {
  VolumeGradientBase::registerKeywords( keys );
  if( keys.reserved("VMEAN") ) keys.use("VMEAN");
  keys.use("MEAN"); keys.use("LESS_THAN"); keys.use("MORE_THAN");
  keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("SUM");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.addFlag("OUTSIDE",false,"calculate quantities for colvars that are on atoms outside the region of interest");
}

ActionVolume::ActionVolume(const ActionOptions&ao):
  Action(ao),
  VolumeGradientBase(ao)
{
  // Find number of quantities
  if( getPntrToMultiColvar()->isDensity() ) nquantities=2;                           // Value + weight
  else if( getPntrToMultiColvar()->getNumberOfQuantities()==2 ) nquantities=2;       // Value + weight
  else nquantities = 1 + getPntrToMultiColvar()->getNumberOfQuantities()-2 + 1;      // Norm  + vector + weight

  // Output some nice information
  std::string functype=getPntrToMultiColvar()->getName();
  std::transform( functype.begin(), functype.end(), functype.begin(), tolower );
  log.printf("  calculating %s inside region of insterest\n",functype.c_str() );

  parseFlag("OUTSIDE",not_in); sigma=0.0;
  if( keywords.exists("SIGMA") ) parse("SIGMA",sigma);
  if( keywords.exists("KERNEL") ) parse("KERNEL",kerneltype);

  if( getPntrToMultiColvar()->isDensity() ) {
    std::string input;
    addVessel( "SUM", input, -1 );  // -1 here means that this value will be named getLabel()
  }
  readVesselKeywords();
}

void ActionVolume::calculateAllVolumes( const unsigned& curr, MultiValue& outvals ) const {
  Vector catom_pos=getPntrToMultiColvar()->getCentralAtomPos( curr );

  double weight; Vector wdf; Tensor vir; std::vector<Vector> refders( getNumberOfAtoms() );
  weight=calculateNumberInside( catom_pos, wdf, vir, refders );
  if( not_in ) {
    weight = 1.0 - weight; wdf *= -1.; vir *=-1;
    for(unsigned i=0; i<refders.size(); ++i) refders[i]*=-1;
  }
  setNumberInVolume( 0, curr, weight, wdf, vir, refders, outvals );
}

bool ActionVolume::inVolumeOfInterest( const unsigned& curr ) const {
  Vector catom_pos=getPntrToMultiColvar()->getCentralAtomPos( curr );
  Vector wdf; Tensor vir; std::vector<Vector> refders( getNumberOfAtoms() );
  double weight=calculateNumberInside( catom_pos, wdf, vir, refders );
  if( not_in ) weight = 1.0 - weight;
  if( weight<getTolerance() ) return false;
  return true;
}

}
}
