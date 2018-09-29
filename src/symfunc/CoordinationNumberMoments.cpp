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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR COORDINATION_MOMENTS
/*

*/
//+ENDPLUMEDOC


class CoordinationNumberMoments : public SymmetryFunctionBase {
private:
  unsigned r_power;
public:
  static void registerKeywords( Keywords& keys );
  explicit CoordinationNumberMoments(const ActionOptions&);
  void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(CoordinationNumberMoments,"COORDINATION_MOMENTS_MATINP")

void CoordinationNumberMoments::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys );
  keys.add("compulsory","R_POWER","the power to which you want to raise the distance");
}

CoordinationNumberMoments::CoordinationNumberMoments(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  parse("R_POWER",r_power);
  log.printf("  multiplying switching function by r^%d\n", r_power);
  addValueWithDerivatives();
}

void CoordinationNumberMoments::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  double d = distance.modulo(); double raised = pow( d, r_power - 1 );
  addToValue( 0, val*raised*d, myvals );
  addWeightDerivative( 0, d*raised, myvals );
  addVectorDerivatives( 0, val*r_power*distance, myvals );
}

class CoordinationNumberMomentsShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit CoordinationNumberMomentsShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CoordinationNumberMomentsShortcut,"COORDINATION_MOMENTS")

void CoordinationNumberMomentsShortcut::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::shortcutKeywords( keys );
}

CoordinationNumberMomentsShortcut::CoordinationNumberMomentsShortcut(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  SymmetryFunctionBase::createSymmetryFunctionObject( getShortcutLabel(), "COORDINATION_MOMENTS", false, false, this );
} 

}
}
