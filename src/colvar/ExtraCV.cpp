/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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
#include "Colvar.h"
#include "ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/Atoms.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR EXTRACV
/*
XX

\par Examples

XX

*/
//+ENDPLUMEDOC


class ExtraCV : public Colvar {
  std::string name;
public:
  explicit ExtraCV(const ActionOptions&);
// active methods:
  void prepare() override;
  void calculate() override;
  unsigned getNumberOfDerivatives() override;
  static void registerKeywords( Keywords& keys );
};


using namespace std;


PLUMED_REGISTER_ACTION(ExtraCV,"EXTRACV")

ExtraCV::ExtraCV(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  addValueWithDerivatives(); setNotPeriodic();
  getPntrToValue()->resizeDerivatives(1);
  parse("NAME",name);
  log<<"  name: "<<name<<"\n";
  isExtraCV=true;
  setExtraCV(name);
}

void ExtraCV::registerKeywords( Keywords& keys ) {
  Action::registerKeywords( keys );
  ActionAtomistic::registerKeywords( keys );
  ActionWithValue::registerKeywords( keys );
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("compulsory","NAME","name of the CV as computed by the MD engine");
}

unsigned ExtraCV::getNumberOfDerivatives() {
  return 1;
}

void ExtraCV::prepare() {
/// \todo: notify Atoms that this is requested
}

// calculator
void ExtraCV::calculate() {
  double value=plumed.getAtoms().getExtraCV(name);
  setValue( value );
  getPntrToComponent(0)->addDerivative(0,1.0);
}

}
}



