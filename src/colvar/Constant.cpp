/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2016 The plumed team
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

namespace PLMD{
namespace colvar{

//+PLUMEDOC COLVAR CONSTANT
/*
Return one or more constant quantities
with or without derivatives.

Useful in combination with functions that
takes in input constants or parameters.

\par Examples

The following input instructs plumed to compute the distance
between atoms 1 and 2. If this distance is between 1.0 and 2.0, it is
printed. If it is lower than 1.0 (larger than 2.0), 1.0 (2.0) is printed

\verbatim
cn: CONSTANT VALUES=1.0,2.0
dis: DISTANCE ATOMS=1,2
sss: SORT ARG=cn.v_0,dis,cn.v_1
PRINT ARG=sss.2
\endverbatim
(See also \ref DISTANCE, \ref SORT, and \ref PRINT).

In case you want to pass a single value you can use VALUE:
\verbatim
cn: CONSTANT VALUE=1.0
dis: DISTANCE ATOMS=1
sss: SORT ARG=cn,dis
PRINT ARG=sss.1
\endverbatim

*/
//+ENDPLUMEDOC

using namespace std;

class Constant : public Colvar {
  vector<double> values;
public:
  explicit Constant(const ActionOptions&);
  virtual void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Constant,"CONSTANT")

Constant::Constant(const ActionOptions&ao):
PLUMED_COLVAR_INIT(ao)
{
  bool noderiv=false;
  parseFlag("NODERIV",noderiv);
  parseVector("VALUES",values);
  if(values.size()==0){
    double v;
    parse("VALUE",v);
// this checks if v is different from NAN
    if(v*2!=v || v==0.0){
      values.resize(1);
      values[0]=v;
    }
  }
  if(values.size()==0) error("Either VALUE or VALUES should be used with CONSTANT");
  checkRead();
  if(values.size()==1) {
    if(!noderiv) addValueWithDerivatives();
    else addValue();
    setNotPeriodic();
    setValue(values[0]);
  } else if(values.size()>1) {
    for(unsigned i=0;i<values.size();i++) {
      std::string num; Tools::convert(i,num);
      if(!noderiv) addComponentWithDerivatives("v_"+num);
      else addComponent("v_"+num);
      componentIsNotPeriodic("v_"+num);
      Value* comp=getPntrToComponent("v_"+num);
      comp->set(values[i]);
    }
  }
// fake request to avoid errors:
  std::vector<AtomNumber> atoms;
  requestAtoms(atoms);
}

void Constant::registerKeywords( Keywords& keys ){
  Colvar::registerKeywords( keys );
  componentsAreNotOptional(keys);
  useCustomisableComponents(keys);
  keys.remove("NUMERICAL_DERIVATIVES"); 
  keys.add("compulsory","VALUES","NAN","The values of the constants");
  keys.add("compulsory","VALUE","NAN","The value of the constant");
  keys.addFlag("NODERIV",false,"Set to TRUE if you want values without derivatives.");  
  keys.addOutputComponent("v","default","the # value"); 
}

// calculator
void Constant::calculate(){
  if(values.size()==1) { 
    setValue(values[0]);
    return;
  }
  for(unsigned i=0;i<values.size();i++) {
    Value* comp=getPntrToComponent(i);
    comp->set(values[i]);
  }
}

}
}



