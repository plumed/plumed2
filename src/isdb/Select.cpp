/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
/* 

*/
#include "function/Function.h"
#include "function/ActionRegister.h"
#include "core/PlumedMain.h"
#include <string>

using namespace std;

namespace PLMD{
namespace isdb{

//+PLUMEDOC FUNCTION SELECT
/*

*/
//+ENDPLUMEDOC

class Select : public function::Function
{
 string selector_;

public:
  explicit Select(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(Select,"SELECT")

void Select::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","SELECTOR","name of the variable used to select");
}

Select::Select(const ActionOptions&ao):
Action(ao), Function(ao)
{
  // name of selector 
  parse("SELECTOR", selector_);

  addValueWithDerivatives(); setNotPeriodic();
  checkRead();

  log.printf("  select based on %s\n",selector_.c_str()); 
}

void Select::calculate()
{
  unsigned iselect = static_cast<unsigned>(plumed.passMap[selector_]); 

  // check if iselect is smaller than the number of arguments
  if(iselect>=getNumberOfArguments()) error("the value of the SELECTOR is greater than the number of arguments!");
  
  // put all the derivatives to zero 
  for(unsigned i=0; i<getNumberOfArguments(); ++i) setDerivative(i, 0.0);

  // set value and derivative for selected argument
  setValue(getArgument(iselect));
  setDerivative(iselect, 1.0);
}

}
}

