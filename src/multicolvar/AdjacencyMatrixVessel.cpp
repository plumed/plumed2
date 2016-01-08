/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2015 The plumed team
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
#include "vesselbase/ActionWithVessel.h"
#include "AdjacencyMatrixVessel.h"
#include "AdjacencyMatrixAction.h" 

namespace PLMD {
namespace multicolvar {

void AdjacencyMatrixVessel::registerKeywords( Keywords& keys ){
  StoreDataVessel::registerKeywords(keys);
}

AdjacencyMatrixVessel::AdjacencyMatrixVessel( const vesselbase::VesselOptions& da ):
StoreDataVessel(da)
{
  function=dynamic_cast<AdjacencyMatrixAction*>( getAction() );
  plumed_assert( function );
}

void AdjacencyMatrixVessel::prepare(){
  finished=false; 
  StoreDataVessel::prepare();
}

void AdjacencyMatrixVessel::setFinishedTrue(){
  finished=true;
}

bool AdjacencyMatrixVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  if( !finished ) return StoreDataVessel::calculate( current, myvals, buffer, der_list );
  return false;
}

void AdjacencyMatrixVessel::finish( const std::vector<double>& buffer ){
  if( !finished ){
     finished=true;
     StoreDataVessel::finish( buffer ); 
     function->dertime=true;
     function->completeCalculation();
  }
}

}
}

