/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "Sum.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "core/ActionRegister.h"

namespace PLMD {
namespace function {

typedef FunctionShortcut<Sum> SumShortcut;
PLUMED_REGISTER_ACTION(SumShortcut,"SUM")
PLUMED_REGISTER_ACTION(SumShortcut,"MEAN")
typedef FunctionOfScalar<Sum> ScalarSum;
PLUMED_REGISTER_ACTION(ScalarSum,"SUM_SCALAR")
PLUMED_REGISTER_ACTION(ScalarSum,"MEAN_SCALAR")
typedef FunctionOfVector<Sum> VectorSum;
PLUMED_REGISTER_ACTION(VectorSum,"SUM_VECTOR")
PLUMED_REGISTER_ACTION(VectorSum,"MEAN_VECTOR")

void Sum::registerKeywords( Keywords& keys ) {
  keys.use("PERIODIC");
}
 
void Sum::read( ActionWithArguments* action ) {
  if( action->getNumberOfArguments()!=1 ) action->error("should only be one argument to sum actions");
}

void Sum::setPrefactor( ActionWithArguments* action, const double pref ) {
  if(action->getName().find("MEAN")!=std::string::npos) prefactor = pref / (action->getPntrToArgument(0))->getNumberOfValues();
  else prefactor = pref;
}

bool Sum::zeroRank() const {
  return true;
}

void Sum::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  vals[0]=prefactor*args[0]; derivatives(0,0)=prefactor;
}

}
}
