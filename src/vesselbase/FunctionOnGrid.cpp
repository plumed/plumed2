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
#include "FunctionOnGrid.h"
#include "VesselRegister.h"

namespace PLMD{
namespace vesselbase{

PLUMED_REGISTER_VESSEL(FunctionOnGrid,"GRID_NOSPLINE")

void FunctionOnGrid::reserveKeyword( Keywords& keys ){
  keys.reserve("optional","GRID_NOSPLINE","create a grid to store a function");
} 

void FunctionOnGrid::registerKeywords( Keywords& keys ){
  GridVesselBase::registerKeywords( keys );
}

FunctionOnGrid::FunctionOnGrid( const VesselOptions& da ):
GridVesselBase(da)
{
  std::string num;
  std::vector<std::string> names(dimension+1);  
  for(unsigned i=0;i<dimension;++i){ Tools::convert(i+1,num); names[i]="x" + num; }  
  names[dimension]=getAction()->getLabel();
  std::vector<bool> mypbc( dimension, false );
  finishSetup( 1, mypbc, names );
}

std::string FunctionOnGrid::description(){
  return getGridDescription();
}

bool FunctionOnGrid::calculate(){
  plumed_merror("This should not be called");
  return true;
}

void FunctionOnGrid::finish(){
  plumed_merror("This should not be called");
}

bool FunctionOnGrid::applyForce( std::vector<double>& forces ){
  plumed_merror("This should not be called");
  return false;
}

}
}
