/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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
#include "core/ActionRegister.h"
#include <string>
#include <algorithm>

using namespace std;

namespace PLMD{
namespace generic{

//+PLUMEDOC FUNCTION PARAMETER
/*
Transform an input argument into a constant value.

\par Examples
The following input tells plumed to calculate a weighted distance where 
the weight is derived from a bias

\verbatim
d: DISTANCE LABEL=dist ATOMS=3,5
res: RESTRAINT ARG=d AT=1 KAPPA=1
p: PARAMETER ARG=res.bias
w: MATHEVAL ARG=d,p.res_bias_par VAR=x,a FUNC=d*exp(+a/2.49) PERIODIC=NO
PRINT ARG=d,w
\endverbatim
(See also \ref DISTANCE, \ref RESTRAINT, \ref MATHEVAL, \ref PRINT).


*/
//+ENDPLUMEDOC


class Parameter :
  public ActionWithValue,
  public ActionWithArguments
{
public:
  Parameter(const ActionOptions&);
  virtual void calculate();
  virtual void apply(){};
  static void registerKeywords(Keywords& keys);
  unsigned getNumberOfDerivatives(){ return 0; }
};


PLUMED_REGISTER_ACTION(Parameter,"PARAMETER")

void Parameter::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords( keys );
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG"); 
  componentsAreNotOptional(keys);
  keys.addOutputComponent("_par","default","one or multiple instances of this quantity will be refereceable elsewhere in the input file. "
                                            "these quantities will named with  the arguments of the bias followed by "
                                            "the character string _bias. These quantities tell the user how much the bias is "
                                            "due to each of the colvars.");

}

Parameter::Parameter(const ActionOptions&ao):
Action(ao),
ActionWithValue(ao),
ActionWithArguments(ao)
{
  checkRead();
  for(unsigned i=0;i<getNumberOfArguments();++i){
     string ss=getPntrToArgument(i)->getName()+"_par";
     addComponentWithDerivatives(ss); 
     componentIsNotPeriodic(ss);
     getPntrToComponent(i)->resizeDerivatives(1);
  }
  clearDependencies(); 
}

void Parameter::calculate(){
  for(unsigned i=0;i<getNumberOfArguments();++i){
     Value* v=getPntrToComponent(i);
     double cv = getArgument(i);
     v->set(cv);
     //setDerivative(v,i,0);
  }
}

}
}


