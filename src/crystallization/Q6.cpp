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
#include "Steinhardt.h"
#include "LocalSteinhardt.h"
#include "core/ActionRegister.h"

//+PLUMEDOC MCOLVAR Q6
/*
Calculate 6th order Steinhardt parameters.

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVARF LOCAL_Q6
/*
Calculate 4th order Steinhardt parameters.

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class Q6 : public Steinhardt {
public:
  static void registerKeywords( Keywords& keys );
  Q6( const ActionOptions& ao );
};

PLUMED_REGISTER_ACTION(Q6,"Q6")
typedef LocalSteinhardt<Q6> LOCAL_Q6;
PLUMED_REGISTER_ACTION(LOCAL_Q6,"LOCAL_Q6")

void Q6::registerKeywords( Keywords& keys ){
  Steinhardt::registerKeywords( keys );
}

Q6::Q6(const ActionOptions& ao ):
Action(ao),
Steinhardt(ao)
{
  setAngularMomentum(6);

  normaliz.resize( 7 );
  normaliz[0] = sqrt( ( 13.0*720.0 ) / (4.0*pi*720.0) );
  normaliz[1] = -sqrt( ( 13.0*120.0 ) / (4.0*pi*5040) );
  normaliz[2] = sqrt( ( 13.0*24) / (4.0*pi*40320) );
  normaliz[3] = -sqrt( ( 13.0*6) / (4.0*pi*362880) );
  normaliz[4] = sqrt( (13.0*2) / (4.0*pi*3628800) );
  normaliz[5] = -sqrt( (13.0*1) / (4.0*pi*39916800) );
  normaliz[6] = sqrt( (13.0*1) / (4.0*pi*479001600) );

  coeff_poly.resize( 7 ); 
  coeff_poly[0]=-0.3125; coeff_poly[1]=0.0;
  coeff_poly[2]=6.5625; coeff_poly[3]=0.0;
  coeff_poly[4]=-19.6875; coeff_poly[5]=0.0;
  coeff_poly[6]=14.4375; 
}

}
}

