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
#include "Colvar.h"
#include "ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR TEMPLATE
/*
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->

\par Examples

<!---You should put an example of how to use your CV here--->

\plumedfile
# This should be a sample input.
t: TEMPLATE ATOMS=1,2
PRINT ARG=t STRIDE=100 FILE=COLVAR
\endplumedfile
<!---You should reference here the other actions used in this example--->
(see also \ref PRINT)

*/
//+ENDPLUMEDOC

class Template : public Colvar {
  bool pbc;

public:
  explicit Template(const ActionOptions&);
// active methods:
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Template,"TEMPLATE")

void Template::registerKeywords(Keywords& keys) {
  Colvar::registerKeywords(keys);
  keys.addFlag("TEMPLATE_DEFAULT_OFF_FLAG",false,"flags that are by default not performed should be specified like this");
  keys.addFlag("TEMPLATE_DEFAULT_ON_FLAG",true,"flags that are by default performed should be specified like this");
  keys.add("compulsory","TEMPLATE_COMPULSORY","all compulsory keywords should be added like this with a description here");
  keys.add("optional","TEMPLATE_OPTIONAL","all optional keywords that have input should be added like a description here");
  keys.add("atoms","ATOMS","the keyword with which you specify what atoms to use should be added like this");
}

Template::Template(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=2)
    error("Number of specified atoms should be 2");
  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;
  checkRead();

  log.printf("  between atoms %d %d\n",atoms[0].serial(),atoms[1].serial());
  if(pbc) log.printf("  using periodic boundary conditions\n");
  else    log.printf("  without periodic boundary conditions\n");

  addValueWithDerivatives(); setNotPeriodic();

  requestAtoms(atoms);
}


// calculator
void Template::calculate() {

  Vector distance;
  if(pbc) {
    distance=pbcDistance(getPosition(0),getPosition(1));
  } else {
    distance=delta(getPosition(0),getPosition(1));
  }
  const double value=distance.modulo();
  const double invvalue=1.0/value;

  setAtomsDerivatives(0,-invvalue*distance);
  setAtomsDerivatives(1,invvalue*distance);
  setBoxDerivatives  (-invvalue*Tensor(distance,distance));
  setValue           (value);
}

}
}



