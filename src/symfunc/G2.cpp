/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "SymmetryFunctionBase.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionRegister.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR G2
/*

\par Examples

*/
//+ENDPLUMEDOC


class G2 : public SymmetryFunctionBase {
private:
  double center, nu;
public:
  static void registerKeywords( Keywords& keys );
  explicit G2(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(G2,"G2")

void G2::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys ); 
  keys.add("compulsory","CENTER","position of the gaussian center");
  keys.add("compulsory","NU","value of the width parameter for the gaussian");
}

G2::G2(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  parse("NU",nu); parse("CENTER",center);
  log.printf("  gaussian center is as %f and width paramter is %f \n", center, nu );
  addValueWithDerivatives();
}

void G2::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  double dlen = distance.modulo(), diff = dlen - center, ee = exp( - nu*diff*diff );
  addToValue( 0, ee*val, myvals ); addWeightDerivative( 0, ee, myvals ); 
  addVectorDerivatives( 0, -(val/dlen)*2*nu*diff*ee*distance, myvals );
}

}
}
