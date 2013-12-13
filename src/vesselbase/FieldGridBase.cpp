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
#include "FieldGridBase.h"

namespace PLMD{
namespace vesselbase{

void FieldGridBase::registerKeywords( Keywords& keys ){
  GridVesselBase::registerKeywords( keys );
}

FieldGridBase::FieldGridBase( const VesselOptions& da ):
GridVesselBase(da)
{
  derlow.resize( dimension );
}

void FieldGridBase::resize(){
  GridVesselBase::resize();
  forces.resize( getAction()->getNumberOfDerivatives() );
}

std::string FieldGridBase::getBaseCVName( const unsigned& icv ) const {
  return getQuantityDescription( dimension+(icv+1)*(1+dimension) );
}

void FieldGridBase::accumulate( const double& val, const double& lowdv, const double& highdv, const double& crossder, const unsigned& derno ){
  // Retrieve the derivative in the low dimenisonal space
  for(unsigned i=0;i<dimension;++i) derlow[i]=getAction()->getElementDerivative(i);
  // Add the value to the grid value
  addToGridElement( currentGridPoint, 0, val );
  // Add the derivatives in the low-dimenisonal space to the grid
  for(unsigned i=0;i<dimension;++i) addToGridElement( currentGridPoint, i+1, lowdv*derlow[i] );
  // Now merge the derivatives 
  // Merge the derivative of the stress
  getAction()->chainRuleForElementDerivatives( currentGridPoint, derno, dimension+1, 0, highdv, this );
  // Merge the derivatives of the stress wrt each property
  for(unsigned i=0;i<dimension;++i){
      getAction()->chainRuleForElementDerivatives( currentGridPoint, derno, dimension+1, i+1, crossder*derlow[i], this );
  }
}

bool FieldGridBase::applyForce(std::vector<double>& ff){
  plumed_dbg_assert( ff.size()==getAction()->getNumberOfDerivatives() );
  if( wasforced ){
      for(unsigned i=0;i<ff.size();++i) ff[i]=forces[i];
  }
  return wasforced;
}



}
}
