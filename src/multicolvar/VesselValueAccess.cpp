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
#include "ActionWithDistribution.h"

namespace PLMD {

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

VesselStoreAllValues::VesselStoreAllValues( const VesselOptions& da ):
VesselValueAccess(da)
{
  setNumberOfValues( getAction()->getNumberOfFunctionsInAction() );
}

void VesselStoreAllValues::resize(){
  ActionWithDistribution* aa=getAction();
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

VesselAccumulator::VesselAccumulator( const VesselOptions& da ):
VesselValueAccess(da),
nbuffers(0)
{
}

void VesselAccumulator::addBufferedValue(){
   nbuffers++;
}

void VesselAccumulator::addOutput( const std::string& label ){
   ActionWithValue* a=dynamic_cast<ActionWithValue*>( getAction() );
   plumed_massert(a,"cannot create passable values as base action does not inherit from ActionWithValue");

   a->addComponentWithDerivatives( label ); 
   a->componentIsNotPeriodic( label );
   final_values.push_back( a->copyOutput( a->getNumberOfComponents()-1 ) );
   setNumberOfValues( nbuffers + final_values.size() );
}


void VesselAccumulator::resize(){
  unsigned nder=getAction()->getNumberOfDerivatives();
  unsigned nfunc=final_values.size(); std::vector<unsigned> sizes( nbuffers + nfunc );
  for(unsigned i=0;i<nbuffers;++i){ sizes[i]=nder; }
  for(unsigned i=0;i<nfunc;++i){ sizes[nbuffers + i]=nder; final_values[i]->resizeDerivatives( nder ); }
  setValueSizes( sizes );
}

bool VesselAccumulator::applyForce( std::vector<double>& forces ){
  std::vector<double> tmpforce( forces.size() );
  forces.assign(forces.size(),0.0); bool wasforced=false;
  for(unsigned i=0;i<final_values.size();++i){
     if( final_values[i]->applyForce( tmpforce ) ){
         wasforced=true;
         for(unsigned j=0;j<forces.size();++j) forces[j]+=tmpforce[j];
     }
  }
  return wasforced;
}

}
