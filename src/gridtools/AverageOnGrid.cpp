/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2015 The plumed team
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
#include "AverageOnGrid.h"

namespace PLMD {
namespace gridtools {

void AverageOnGrid::registerKeywords( Keywords& keys ){
  HistogramOnGrid::registerKeywords( keys );
}

AverageOnGrid::AverageOnGrid( const vesselbase::VesselOptions& da ):
HistogramOnGrid(da)
{
  arg_names.push_back( "density" );
  if( !discrete ){
    for(unsigned i=0;i<dimension;++i) arg_names.push_back( "ddensity_" + arg_names[i] );
    nper += (dimension+1);
  } else {
    nper += 1;
  }
}

void AverageOnGrid::accumulate( const unsigned& ipoint, const double& weight, const double& dens, const std::vector<double>& der, std::vector<double>& buffer ) const {
  buffer[bufstart+nper*ipoint] += weight*dens; buffer[ bufstart+nper*(ipoint+1) - (dimension+1) ] += dens;
  if( der.size()>0 ){
    for(unsigned j=0;j<dimension;++j) buffer[ bufstart+nper*ipoint + 1 + j ] += weight*der[j];
    for(unsigned j=0;j<dimension;++j) buffer[ bufstart+nper*(ipoint+1) - dimension + j ] += der[j];
  }
}

double AverageOnGrid::getGridElement( const unsigned& ipoint, const unsigned& jelement ) const {
  if( unormalised ) return data[nper*ipoint + jelement];  

  if( noderiv ) return data[nper*ipoint+jelement] / data[nper*(1+ipoint) - 1];

  double rdenom = 1.0; 
  if( fabs(data[nper*(ipoint+1) -(dimension+1)])>epsilon ) rdenom = 1. / data[nper*(ipoint+1) - (dimension+1)];
  
  unsigned jderiv = jelement%(1+dimension);
  if( jderiv==0 ) return rdenom*data[nper*ipoint+jelement]; 

  unsigned jfloor = std::floor( jelement / (1+dimension) );
  return rdenom*data[nper*ipoint+jelement] - rdenom*rdenom*data[nper*ipoint+jfloor]*data[nper*(ipoint+1) - (dimension+1) + jderiv ]; 
}

double AverageOnGrid::getGridElementForPrint( const unsigned& ipoint, const unsigned& jelement ) const {
  plumed_assert( bounds_set && ipoint<npoints && jelement<nper && active[ipoint] );
  return data[nper*ipoint + jelement]; 
}

}
}
