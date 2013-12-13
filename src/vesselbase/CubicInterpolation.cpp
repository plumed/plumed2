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
#include "CubicInterpolation.h"

namespace PLMD {
namespace vesselbase {

CubicInterpolation::CubicInterpolation( GridVesselBase* gg, const unsigned dstart ):
InterpolationBase(gg,dstart)
{
  plumed_dbg_assert( getDimension()==1 );
  clist.resize( 4*getNumberOfSplinePoints() );
} 

void CubicInterpolation::setInterpolationTables(){
  std::vector<unsigned> np(1); std::vector<double> der(1);
  getGridPointSpacing(der); double norm = (der[0]*der[0])/6.0;
  for(unsigned i=0;i<getNumberOfSplinePoints()-1;++i){
      unsigned pij=i*4; np[0]=i; 
      clist[pij] = getValueAndDerivatives( np, der ); 
      clist[pij+2] = der[0]*norm; np[0]=i+1; 
      clist[pij+1] = getValueAndDerivatives( np, der ); 
      clist[pij+3] = der[0]*norm;
  }
}

double CubicInterpolation::interpolateFunction( const unsigned& mybox, const std::vector<double>& dd ){
  plumed_dbg_assert( dd.size()==1 );
  double b = dd[0], a = 1 - b;
  double *cbase=&clist[(mybox*4)+3], *c3=cbase-1, *c2=c3-1, *c1=c2-1;
  double f=a*(*c1) + b*(*c2) + (a*a*a-a)*(*c3) + (b*b*b-b)*(*cbase);
  delete cbase; delete c3; delete c2; delete c1;
  return f;
}

}
}
