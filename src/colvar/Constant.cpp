/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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

//+PLUMEDOC COLVAR CONSTANT
/*
Return one or more constant quantities with or without derivatives.

Useful in combination with functions that
takes in input constants or parameters.

\par Examples

The following input instructs plumed to compute the distance
between atoms 1 and 2. If this distance is between 1.0 and 2.0, it is
printed. If it is lower than 1.0 (larger than 2.0), 1.0 (2.0) is printed

\plumedfile
cn: CONSTANT VALUES=1.0,2.0
dis: DISTANCE ATOMS=1,2
sss: SORT ARG=cn.v-0,dis,cn.v-1
PRINT ARG=sss.2
\endplumedfile

In case you want to pass a single value you can use VALUE:
\plumedfile
cn: CONSTANT VALUE=1.0
dis: DISTANCE ATOMS=1,2
sss: SORT ARG=cn,dis
PRINT ARG=sss.1
\endplumedfile

*/
//+ENDPLUMEDOC

using namespace std;

class Constant : public Colvar {
  vector<double> values;
public:
  explicit Constant(const ActionOptions&);
  void calculate() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Constant,"CONSTANT")

Constant::Constant(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao)
{
  bool noderiv=false;
  parseFlag("NODERIV",noderiv);
  parseVector("VALUES",values);
  vector<double> value;
  parseVector("VALUE",value);
  if(values.size()==0&&value.size()==0) error("One should use either VALUE or VALUES");
  if(values.size()!=0&&value.size()!=0) error("One should use either VALUE or VALUES");
  if(value.size()>1) error("VALUE cannot take more than one number");
  if(values.size()==0) {
    values.resize(1);
    values[0]=value[0];
  }
  checkRead();
  if(values.size()==1) {
    if(!noderiv) addValueWithDerivatives();
    else addValue();
    setNotPeriodic();
    setValue(values[0]);
  } else if(values.size()>1) {
    for(unsigned i=0; i<values.size(); i++) {
      std::string num; Tools::convert(i,num);
      if(!noderiv) addComponentWithDerivatives("v-"+num);
      else addComponent("v-"+num);
      componentIsNotPeriodic("v-"+num);
      Value* comp=getPntrToComponent("v-"+num);
      comp->set(values[i]);
    }
  }
// fake request to avoid errors:
  std::vector<AtomNumber> atoms;
  requestAtoms(atoms);
}

void Constant::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  componentsAreNotOptional(keys);
  keys.remove("NUMERICAL_DERIVATIVES");
  keys.add("optional","VALUES","The values of the constants");
  keys.add("optional","VALUE","The value of the constant");
  keys.addFlag("NODERIV",false,"Set to TRUE if you want values without derivatives.");
  keys.addOutputComponent("v","default","the # value");
}

// calculator
void Constant::calculate() {
  if(values.size()==1) {
    setValue(values[0]);
    return;
  }
  for(unsigned i=0; i<values.size(); i++) {
    Value* comp=getPntrToComponent(i);
    comp->set(values[i]);
  }
}

}
}



