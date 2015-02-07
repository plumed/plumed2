/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
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
#include "ActionVolume.h"

namespace PLMD {
namespace multicolvar {

void ActionVolume::registerKeywords( Keywords& keys ){
  VolumeGradientBase::registerKeywords( keys );
  if( keys.reserved("VMEAN") ) keys.use("VMEAN");
  keys.use("MEAN"); keys.use("LESS_THAN"); keys.use("MORE_THAN");
  keys.use("BETWEEN"); keys.use("HISTOGRAM");
  keys.add("compulsory","SIGMA","the width of the function to be used for kernel density estimation");
  keys.add("compulsory","KERNEL","gaussian","the type of kernel function to be used");
  keys.addFlag("OUTSIDE",false,"calculate quantities for colvars that are on atoms outside the region of interest");
}

ActionVolume::ActionVolume(const ActionOptions&ao):
Action(ao),
VolumeGradientBase(ao)
{
  // Find number of quantities
  if( getPntrToMultiColvar()->isDensity() ) nquantities=5;                           // Value + catom + weight 
  else if( getPntrToMultiColvar()->getNumberOfQuantities()==5 ) nquantities=5;       // Value + catom + weight
  else nquantities = 1 + 3 + getPntrToMultiColvar()->getNumberOfQuantities()-5 + 1;  // Norm + catom + vector + weight 

  // Output some nice information
  std::string functype=getPntrToMultiColvar()->getName();
  std::transform( functype.begin(), functype.end(), functype.begin(), tolower );
  log.printf("  calculating %s inside region of insterest\n",functype.c_str() ); 

  parseFlag("OUTSIDE",not_in); parse("SIGMA",sigma); 
  bead.isNotPeriodic(); 
  std::string kerneltype; parse("KERNEL",kerneltype); 
  bead.setKernelType( kerneltype );
  
  if( getPntrToMultiColvar()->isDensity() ){
     std::string input;
     addVessel( "SUM", input, -1 );  // -1 here means that this value will be named getLabel()
  } else {
     readVesselKeywords();
  }
}

void ActionVolume::calculateAllVolumes(){
  Vector catom_pos=getPntrToMultiColvar()->retrieveCentralAtomPos();

  double weight; Vector wdf; 
  weight=calculateNumberInside( catom_pos, bead, wdf ); 
  if( not_in ){ weight = 1.0 - weight; wdf *= -1.; }  

  setNumberInVolume( nquantities-1, weight, wdf );
}

}
}
