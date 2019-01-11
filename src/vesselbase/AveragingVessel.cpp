/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2019 The plumed team
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
#include "AveragingVessel.h"
#include "ActionWithAveraging.h"

namespace PLMD {
namespace vesselbase {

void AveragingVessel::registerKeywords( Keywords& keys ) {
  Vessel::registerKeywords( keys );
}

AveragingVessel::AveragingVessel( const vesselbase::VesselOptions& vo ):
  Vessel(vo),
  wascleared(true)
{
  if( getAction() ) {
    ActionWithAveraging* myav = dynamic_cast<ActionWithAveraging*>( getAction() );
    plumed_assert( myav ); unormalised = myav->ignoreNormalization();
  }
}

void AveragingVessel::finish( const std::vector<double>& buffer ) {
  wascleared=false; for(unsigned i=1; i<data.size(); ++i) data[i]+=buffer[bufstart + i - 1];
}

bool AveragingVessel::wasreset() const {
  return wascleared;
}

void AveragingVessel::clear() {
  plumed_assert( wascleared ); data.assign( data.size(), 0.0 );
}

void AveragingVessel::reset() {
  wascleared=true;
}

void AveragingVessel::setDataSize( const unsigned& size ) {
  if( data.size()!=(1+size) ) data.resize( 1+size, 0 );
}

}
}

