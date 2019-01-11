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
#include "AverageVessel.h"

namespace PLMD {
namespace analysis {

void AverageVessel::registerKeywords( Keywords& keys ) {
  vesselbase::AveragingVessel::registerKeywords( keys );
  keys.add("optional","PERIODIC","is the quantity being averaged periodic and what is its domain");
}

AverageVessel::AverageVessel( const vesselbase::VesselOptions& da):
  AveragingVessel(da)
{
  parseVector("PERIODIC",domain);
  plumed_assert( domain.size()==2 || domain.size()==0 );
}

void AverageVessel::resize() {
  resizeBuffer(0);
  if( domain.size()==2 ) setDataSize(2);
  else setDataSize(1);
}

void AverageVessel::accumulate( const double& weight, const double& val ) {
  if( domain.size()==2 ) {
    // Average with Berry Phase
    double tval = 2*pi*( val - domain[0] ) / ( domain[1] - domain[0] );
    addDataElement( 0, weight*sin(tval) ); addDataElement( 1, weight*cos(tval) );
  } else addDataElement( 0, weight*val );
}

double AverageVessel::getAverage() const {
  if( domain.size()==2 ) return domain[0] + (( domain[1] - domain[0] )*atan2( getDataElement(0), getDataElement(1) ) / (2*pi));
  return getDataElement(0);
}

void AverageVessel::calculate( const unsigned& current, MultiValue& myvals, std::vector<double>& buffer, std::vector<unsigned>& der_list ) const {
  plumed_error();
}

}
}
