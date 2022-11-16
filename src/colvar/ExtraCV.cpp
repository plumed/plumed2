/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2020 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR EXTRACV
/*
Allow PLUMED to use collective variables computed in the MD engine.

This feature requires the MD engine to use special instructions to pass to PLUMED the value of
some pre-computed collective variable. Check the documentation of the MD code to find out which
collective variables can be computed and passed to PLUMED. These variables can then be accessed by
name using the EXTRACV action.

\par Examples

This example takes the lambda variable pre-computed in GROMACS and apply to it a restraint to keep
it close to the value 3.
\plumedfile
l: EXTRACV NAME=lambda
RESTRAINT ARG=l KAPPA=10 AT=3
\endplumedfile


*/
//+ENDPLUMEDOC


class ExtraCV :
public ActionWithValue, 
public ActionWithArguments 
{
public:
  explicit ExtraCV(const ActionOptions&);
// active methods:
  void calculate() override;
  void apply();
  unsigned getNumberOfDerivatives() const override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(ExtraCV,"EXTRACV")

ExtraCV::ExtraCV(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao)
{
  std::vector<std::string> argn(1); parse("NAME",argn[0]);
  plumed.cmd("createValue " + argn[0] + ": PUT UNIT=number SHAPE=0 PERIODIC=NO"); 
  std::vector<Value*> arg; ActionWithArguments::interpretArgumentList(argn,plumed.getActionSet(),this,arg); 
  if( arg.size()!=1 || arg[0]->getRank()>0 ) error("input argument is not valid"); 
  log<<"  name: "<<arg[0]->getName()<<"\n"; 
  addValueWithDerivatives(); setNotPeriodic();
  getPntrToValue()->resizeDerivatives(1); requestArguments( arg, false );
}

void ExtraCV::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("compulsory","NAME","name of the CV as computed by the MD engine");
}

unsigned ExtraCV::getNumberOfDerivatives() const {
  return 1;
}

// calculator
void ExtraCV::calculate() {
  setValue( getPntrToArgument(0)->get() );
  getPntrToComponent(0)->addDerivative(0,1.0);
}

void ExtraCV::apply() {
  std::vector<double> forcesToApply(1); unsigned ss=0;
  if( getForcesFromValues( forcesToApply ) ) setForcesOnArguments( 0, forcesToApply, ss );
}

}
}



