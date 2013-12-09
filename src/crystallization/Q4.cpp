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
#include "Steinhardt.h"
#include "LocalSteinhardt.h" 
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR Q4
/*
Calculate 4th order Steinhardt parameters.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVARF LOCAL_Q4
/*
Calculate 4th order Steinhardt parameters.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class Q4 : public Steinhardt {
public:
  static void registerKeywords( Keywords& keys );
  Q4( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(Q4,"Q4")
typedef LocalSteinhardt<Q4> LOCAL_Q4;
PLUMED_REGISTER_ACTION(LOCAL_Q4,"LOCAL_Q4")

void Q4::registerKeywords( Keywords& keys ){
  Steinhardt::registerKeywords( keys );
}

Q4::Q4(const ActionOptions& ao ):
Action(ao),
Steinhardt(ao)
{
  setAngularMomentum(4);

  normaliz.resize( 5 );
  normaliz[0] = sqrt( ( 9.0*24.0 ) / (4.0*pi*24.0) );
  normaliz[1] = -sqrt( ( 9.0*6.0 ) / (4.0*pi*120.0) );
  normaliz[2] = sqrt( ( 9.0*2.0) / (4.0*pi*720.0) );
  normaliz[3] = -sqrt( ( 9.0*1) / (4.0*pi*5040.0) );
  normaliz[4] = sqrt( (9.0*1) / (4.0*pi*40320.0) );

  coeff_poly.resize( 5 ); 
  coeff_poly[0]=0.375; coeff_poly[1]=0.0;
  coeff_poly[2]=-3.75; coeff_poly[3]=0.0;
  coeff_poly[4]=4.375; 
}

}
}

