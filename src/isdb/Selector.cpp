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
#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include <string>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC ISDB_GENERIC SELECTOR
/*
Defines a variable (of the type double) inside the PLUMED code that can be used and modified by other actions.

A \ref SELECTOR can be used for example to activate or modify a bias based on its current value.

\par Examples

A typical example is the simulated-tempering like approach activated by \ref RESCALE.
In this example the total potential energy of the system is scaled
by a parameter defined on a grid of dimension NBIN in the range from 1 to MAX_RESCALE.
The value of the scaling parameter is determined by the current value of the \ref SELECTOR GAMMA.
The value of the \ref SELECTOR is updated by a MC protocol inside the \ref RESCALE class.
A well-tempered metadynamics potential is used to enhance sampling in the \ref SELECTOR space.

\plumedfile
ene:  ENERGY

SELECTOR NAME=GAMMA VALUE=0

RESCALE ...
LABEL=res ARG=ene TEMP=300
SELECTOR=GAMMA MAX_RESCALE=1.2 NBIN=2
W0=1000 BIASFACTOR=100.0 BSTRIDE=2000 BFILE=bias.dat
...

PRINT FILE=COLVAR ARG=* STRIDE=100
\endplumedfile

*/
//+ENDPLUMEDOC

class Selector:
  public Action
{
public:
  static void registerKeywords( Keywords& keys );
  explicit Selector(const ActionOptions&ao);
  void calculate() {}
  void apply() {}
};

PLUMED_REGISTER_ACTION(Selector,"SELECTOR")

void Selector::registerKeywords( Keywords& keys ) {
  Action::registerKeywords(keys);
  keys.add("compulsory","NAME","name of the SELECTOR");
  keys.add("compulsory","VALUE","set (initial) value of the SELECTOR");
}

Selector::Selector(const ActionOptions&ao):
  Action(ao)
{
  string name;
  parse("NAME", name);
  double value;
  parse("VALUE", value);
  plumed.passMap[name] = value;
}

}
}

