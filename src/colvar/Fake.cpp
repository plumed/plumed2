/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2019 The plumed team
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
#include <cassert>
#include <iostream>

using namespace std;

namespace PLMD {
namespace colvar {

//+PLUMEDOC COLVAR FAKE
/*
This is a fake colvar container used by cltools or various other actions that supports input and period definitions

\par Examples

\plumedfile
FAKE ATOMS=1 PERIODIC=-3.14,3.14   LABEL=d2
\endplumedfile

*/
//+ENDPLUMEDOC

class ColvarFake : public Colvar {

public:
  static void registerKeywords( Keywords& keys );
  explicit ColvarFake(const ActionOptions&);
// active methods:
  void calculate() override;
};

PLUMED_REGISTER_ACTION(ColvarFake,"FAKE")

void ColvarFake::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the fake atom index, a number is enough");
  keys.reserve("compulsory","PERIODIC","if the output of your function is periodic then you should specify the periodicity of the function.  If the output is not periodic you must state this using PERIODIC=NO,NO (one for the lower and the other for the upper boundary). For multicomponents then it is PERIODIC=mincomp1,maxcomp1,mincomp2,maxcomp2  etc ");
  keys.use("PERIODIC");
  keys.add("optional","COMPONENTS","additional components that this variable is supposed to have. Periodicity is ruled by PERIODIC keyword ");
  useCustomisableComponents(keys);
}

ColvarFake::ColvarFake(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);

  vector<string> comps;
  // multiple components for this variable
  parseVector("COMPONENTS",comps);
  if(comps.size()!=0) {
    for(unsigned i=0; i<comps.size(); i++) {
      addComponentWithDerivatives(comps[i]);
    }
    // periodicity
  } else {
    // only one component for this variable
    addValueWithDerivatives();
  }
  std::vector<std::string> period;
  parseVector("PERIODIC",period);
  if(period.size()!=0) {
    plumed_massert(static_cast<unsigned>(getNumberOfComponents()*2)==period.size(),"the periodicty should coincide with the number of components");
    if(comps.size()!=0) {
      for(int i=0; i<getNumberOfComponents(); i++) {
        string pp=comps[i];
        if(period[i*2]!="none" && period[i*2+1]!="none" ) {
          componentIsPeriodic(pp,period[i*2],period[i*2+1]);
        } else {
          componentIsNotPeriodic(pp);
        }
      }
    } else {
      if(period[0]!="none" && period[1]!="none" ) {
        setPeriodic(period[0],period[1]);
      } else {
        setNotPeriodic();
      }
    }
  } else {
    if(comps.size()!=0) {
      for(int i=0; i<getNumberOfComponents(); i++) {
        componentIsNotPeriodic(getPntrToComponent(i)->getName());
      }
    } else {
      setNotPeriodic();
    }
  }
  checkRead();
  requestAtoms(atoms);

}


// calculator
void ColvarFake::calculate() {
  plumed_merror("you should never have got here");
}

}
}



