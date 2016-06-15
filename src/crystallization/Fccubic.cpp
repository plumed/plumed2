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
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD{
namespace crystallization{

//+PLUMEDOC MCOLVAR FCCUBIC    
/*

\par Examples

*/
//+ENDPLUMEDOC


class Fccubic : public CubicHarmonicBase {
private:
  double alpha, a1, b1;
public:
  static void registerKeywords( Keywords& keys );
  explicit Fccubic(const ActionOptions&);
  double calculateCubicHarmonic( const Vector& distance, const double& d2, Vector& myder ) const ;
};

PLUMED_REGISTER_ACTION(Fccubic,"FCCUBIC")

void Fccubic::registerKeywords( Keywords& keys ){
  CubicHarmonicBase::registerKeywords( keys );
  keys.add("compulsory","ALPHA","3.0","The alpha parameter of the angular function");
}

Fccubic::Fccubic(const ActionOptions&ao):
Action(ao),
CubicHarmonicBase(ao)
{
  // Scaling factors such that '1' corresponds to fcc lattice
  // and '0' corresponds to isotropic (liquid)
  parse("ALPHA",alpha);
  a1 = 80080. / (2717. + 16*alpha); b1 = 16.*(alpha-143)/(2717+16*alpha);
  log.printf("  setting alpha paramter equal to %f \n",alpha);
  // And setup the ActionWithVessel
  checkRead();
}

double Fccubic::calculateCubicHarmonic( const Vector& distance, const double& d2, Vector& myder ) const {
  double x2 = distance[0]*distance[0];
  double x4 = x2*x2;

  double y2 = distance[1]*distance[1];
  double y4 = y2*y2;

  double z2 = distance[2]*distance[2];
  double z4 = z2*z2;
  
  double r8 = pow( d2, 4 );
  double r12 = pow( d2, 6 );

  double tmp = ((x4*y4)+(x4*z4)+(y4*z4))/r8-alpha*x4*y4*z4/r12;

  double t0 = (x2*y4+x2*z4)/r8-alpha*x2*y4*z4/r12;
  double t1 = (y2*x4+y2*z4)/r8-alpha*y2*x4*z4/r12;
  double t2 = (z2*x4+z2*y4)/r8-alpha*z2*x4*y4/r12;
  double t3 = (2*tmp-alpha*x4*y4*z4/r12)/d2;         
 
  myder[0]=4*a1*distance[0]*(t0-t3);
  myder[1]=4*a1*distance[1]*(t1-t3);
  myder[2]=4*a1*distance[2]*(t2-t3);
         
  return a1*tmp+b1;
}

}
}

