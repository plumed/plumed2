/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013 The plumed team
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
#include "NearestNeighborInterpolation.h"

namespace PLMD {
namespace vesselbase {

NearestNeighborInterpolation::NearestNeighborInterpolation( GridVesselBase* gg, const unsigned dstart ):
InterpolationBase(gg,dstart)
{
  mypows.resize( getDimension() );
  for(unsigned i=0;i<mypows.size();++i) mypows[i]=getGridStride( i );
} 

double NearestNeighborInterpolation::interpolateFunction( const unsigned& mybox, const std::vector<double>& dd ){
  unsigned pp=0;
  for(unsigned i=0;i<dd.size();++i){
      if( dd[i]>0.5 ) pp+=mypows[i];
  }
  return getValue( mybox + pp );
}

}
}
