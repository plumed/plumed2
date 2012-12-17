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
#include "VesselValueAccess.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase{

VesselValueAccess::VesselValueAccess( const VesselOptions& da ) :
Vessel(da)
{
}

void VesselValueAccess::setNumberOfValues( const unsigned& n ){
   value_starts.resize( n + 1 );
}

void VesselValueAccess::setValueSizes( const std::vector<unsigned>& val_sizes ){
  plumed_assert( (val_sizes.size()+1)==value_starts.size() );
  unsigned vstart=0;
  for(unsigned i=0;i<val_sizes.size();++i){ value_starts[i]=vstart; vstart+=val_sizes[i]+1; }
  value_starts[val_sizes.size()]=vstart;
  resizeBuffer( vstart ); 
}

}
}
