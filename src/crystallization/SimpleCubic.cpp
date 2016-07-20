/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2016 The plumed team
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
#include "CubicHarmonicBase.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace crystallization{

//+PLUMEDOC MCOLVAR SIMPLECUBIC
/*
Calculate whether or not the coordination spheres of atoms are arranged as they would be in a simple
cubic structure.

\par Examples

The following input tells plumed to calculate the simple cubic parameter for the atoms 1-100 with themselves.
The mean value is then calculated.
\verbatim
SIMPLECUBIC SPECIES=1-100 R_0=1.0 MEAN
\endverbatim

The following input tells plumed to look at the ways atoms 1-100 are within 3.0 are arranged about atoms
from 101-110.  The number of simple cubic parameters that are greater than 0.8 is then output
\verbatim
SIMPLECUBIC SPECIESA=101-110 SPECIESB=1-100 R_0=3.0 MORE_THAN={RATIONAL R_0=0.8 NN=6 MM=12 D_0=0}
\endverbatim

*/
//+ENDPLUMEDOC


class SimpleCubic : public CubicHarmonicBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit SimpleCubic(const ActionOptions&);
  double calculateCubicHarmonic( const Vector& distance, const double& d2, Vector& myder ) const ;
};

PLUMED_REGISTER_ACTION(SimpleCubic,"SIMPLECUBIC")

void SimpleCubic::registerKeywords( Keywords& keys ){
  CubicHarmonicBase::registerKeywords( keys );
}

SimpleCubic::SimpleCubic(const ActionOptions&ao):
Action(ao),
CubicHarmonicBase(ao)
{
  checkRead();
}

double SimpleCubic::calculateCubicHarmonic( const Vector& distance, const double& d2, Vector& myder ) const { 
  double x2 = distance[0]*distance[0];
  double x3 = distance[0]*x2;
  double x4 = distance[0]*x3;

  double y2 = distance[1]*distance[1];
  double y3 = distance[1]*y2;
  double y4 = distance[1]*y3;         
  
  double z2 = distance[2]*distance[2];
  double z3 = distance[2]*z2; 
  double z4 = distance[2]*z3;

  double r4 = pow( d2, 2 );
  double tmp = ( x4 + y4 + z4 ) / r4;

  double t1=(x2+y2+z2), t2=t1*t1, t3=(x4+y4+z4)/(t1*t2);
  myder[0] = 4*x3/t2-4*distance[0]*t3; 
  myder[1] = 4*y3/t2-4*distance[1]*t3; 
  myder[2] = 4*z3/t2-4*distance[2]*t3; 
  return tmp; 
}

}
}

