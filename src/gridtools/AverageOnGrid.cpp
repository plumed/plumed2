/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2020 The plumed team
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
#include "AverageOnGrid.h"

namespace PLMD {
namespace gridtools {

void AverageOnGrid::registerKeywords( Keywords& keys ) {
  HistogramOnGrid::registerKeywords( keys );
}

AverageOnGrid::AverageOnGrid( const vesselbase::VesselOptions& da ):
  HistogramOnGrid(da)
{
  arg_names.push_back( "density" );
  if( !discrete ) {
    for(unsigned i=0; i<dimension; ++i) arg_names.push_back( "ddensity_" + arg_names[i] );
    nper += (dimension+1);
  } else {
    nper += 1;
  }
}

void AverageOnGrid::accumulate( const unsigned& ipoint, const double& weight, const double& dens, const std::vector<double>& der, std::vector<double>& buffer ) const {
  buffer[bufstart+nper*ipoint] += weight*dens; buffer[ bufstart+nper*(ipoint+1) - (dimension+1) ] += dens;
  if( der.size()>0 ) {
    for(unsigned j=0; j<dimension; ++j) buffer[ bufstart+nper*ipoint + 1 + j ] += weight*der[j];
    for(unsigned j=0; j<dimension; ++j) buffer[ bufstart+nper*(ipoint+1) - dimension + j ] += der[j];
  }
}

double AverageOnGrid::getGridElement( const unsigned& ipoint, const unsigned& jelement ) const {
  if( noAverage() ) return getDataElement( nper*ipoint + jelement);

  if( jelement>=(nper-(dimension+1)) ) return getDataElement( nper*ipoint + jelement );

  if( noderiv ) return getDataElement( nper*ipoint+jelement ) / getDataElement( nper*(1+ipoint) - 1);

  double rdenom = 1.0;
  if( fabs(getDataElement( nper*(ipoint+1) -(dimension+1) ))>epsilon ) rdenom = 1. / getDataElement( nper*(ipoint+1) - (dimension+1) );

  unsigned jderiv = jelement%(1+dimension);
  if( jderiv==0 ) return rdenom*getDataElement( nper*ipoint+jelement );

  unsigned jfloor = std::floor( jelement / (1+dimension) );
  return rdenom*getDataElement( nper*ipoint+jelement ) - rdenom*rdenom*getDataElement(nper*ipoint+jfloor)*getDataElement(nper*(ipoint+1) - (dimension+1) + jderiv);
}

}
}
