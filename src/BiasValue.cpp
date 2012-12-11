/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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
#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS BIASVALUE 
/*
Takes the value of one variable and use it as a bias


\par Examples
The following input tells plumed to restrain the distance between atoms 3 and 5
and the distance between atoms 2 and 4, at different equilibrium
values. It then tells plumed to print the energy of the restraint
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=3,6 LABEL=d2
BIASVALUE ARG=d1,d2 LABEL=b
PRINT ARG=d1,d2,b.d1,b.d2
\endverbatim
(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC

class BiasValue : public Bias{
public:
  BiasValue(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(BiasValue,"BIASVALUE")

void BiasValue::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
}

BiasValue::BiasValue(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao)
{
  checkRead();
  // add one bias for each argument  
 // addComponent("bias");
  for(unsigned i=0;i<getNumberOfArguments();++i){ 
	//log<<getPntrToArgument(i)->getName()<<"\n";
        string ss;
        ss="bias."+getPntrToArgument(i)->getName();
	addComponent(ss); componentIsNotPeriodic(ss);
  }
}

void BiasValue::calculate(){
for(unsigned i=0;i< getNumberOfArguments() ;++i){
  double val; val=getArgument(i); 
//  log<<"BIAS "<<val<<"\n";
  getPntrToComponent(i)->set(val);
  setOutputForce(i,-1.);
}
}

}
