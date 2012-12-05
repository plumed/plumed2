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
#include "ActionWithVessel.h"
#include "VesselStoreAllValues.h"

namespace PLMD {

VesselStoreAllValues::VesselStoreAllValues( const VesselOptions& da ):
VesselValueAccess(da)
{
  setNumberOfValues( getAction()->getNumberOfFunctionsInAction() );
}

void VesselStoreAllValues::resize(){
  ActionWithVessel* aa=getAction();
  unsigned nfunc=aa->getNumberOfFunctionsInAction();
  std::vector<unsigned> sizes( nfunc );
  for(unsigned i=0;i<nfunc;++i) sizes[i]=aa->getNumberOfDerivatives(i);
  setValueSizes( sizes ); local_resizing();
}

bool VesselStoreAllValues::calculate( const unsigned& i, const double& tolerance ){
  Value myvalue=getAction()->retreiveLastCalculatedValue();
  unsigned ider=value_starts[i]; setBufferElement( ider, myvalue.get() ); ider++;
  for(unsigned j=0;j<myvalue.getNumberOfDerivatives();++j){ setBufferElement( ider, myvalue.getDerivative(j) ); ider++; }
  //setValue( i, myvalue );
  return true;
}

}


