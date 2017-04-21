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
#include "core/Action.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include <string>

using namespace std;

namespace PLMD {
namespace isdb {

//+PLUMEDOC GENERIC SELECTOR
/*

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

