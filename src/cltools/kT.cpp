/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
#include "CLTool.h"
#include "CLToolRegister.h"
#include "tools/Tools.h"
#include "config/Config.h"
#include "tools/Units.h"
#include "core/ActionRegister.h"
#include <cstdio>
#include <string>
#include <vector>
#include <iostream>

using namespace std;

namespace PLMD {
namespace cltools {

//+PLUMEDOC TOOLS kt
/*
Print out the value of \f$k_BT\f$ at a particular temperature

\par Examples

The following command will tell you the value of \f$k_BT\f$ when T is equal
to 300 K in eV

\verbatim
plumed kt --temp 300 --units eV
\endverbatim

*/
//+ENDPLUMEDOC

class kt:
  public CLTool
{
public:
  static void registerKeywords( Keywords& keys );
  explicit kt(const CLToolOptions& co );
  int main(FILE* in, FILE*out,Communicator& pc) override;
  string description()const override {
    return "print out the value of kT at a particular temperature";
  }
};

PLUMED_REGISTER_CLTOOL(kt,"kt")

void kt::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--temp","print the manual for this particular action");
  keys.add("compulsory","--units","kj/mol","the units of energy can be kj/mol, kcal/mol, j/mol, eV or the conversion factor from kj/mol");
}

kt::kt(const CLToolOptions& co ):
  CLTool(co)
{
  inputdata=commandline;
}

int kt::main(FILE* in, FILE*out,Communicator& pc) {

  std::string unitname; parse("--units",unitname);
  Units units; units.setEnergy( unitname );
  double temp; parse("--temp",temp);
  double kk=(kBoltzmann*temp)/units.getEnergy();
  std::fprintf(out,"When the temperature is %f kelvin kT is equal to %f %s\n",temp,kk,unitname.c_str());
  return 0;
}

} // End of namespace
}
