/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2015 The plumed team
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
#include "core/ActionPilot.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

using namespace std;

namespace PLMD{
namespace generic{

//+PLUMEDOC GENERIC STOPWHEN
/*
This action can be used in methods like forward flux sampling.  It will terminate your calculation once
the value of a CV or set of CVs is either greater than or less than a certain target value or set of target
values.

\par Examples

The following example will terminate the calculation once the value of the distance between atoms 1 and 2 is
greater than 1 nm

\verbatim
d1: DISTANCE ATOMS=1,2 
STOPWHEN ARG=d1 CONDITION=GT AT=1.0
\endverbatim

*/
//+ENDPLUMEDOC

class StopWhen :
public ActionPilot,
public ActionWithArguments
{
private:
  std::vector<double> targets;
  std::vector<bool> cond;
public:
  void calculate(){}
  explicit StopWhen(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply(){}
  void update();
};

PLUMED_REGISTER_ACTION(StopWhen,"STOPWHEN")

void StopWhen::registerKeywords(Keywords& keys){
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("hidden","STRIDE","1","this should never be changed!!!!");
  keys.add("compulsory","CONDITION","do we stop the calculation when the CV is greater than (GT) or les than (LT) the AT value");
  keys.add("compulsory","AT","the values of the CV at which to stop the calculation");
}

StopWhen::StopWhen(const ActionOptions&ao):
Action(ao),
ActionPilot(ao),
ActionWithArguments(ao),
targets(getNumberOfArguments()),
cond(getNumberOfArguments())
{

  parseVector("AT",targets);
  std::vector<std::string> conditions; parseVector("CONDITION",conditions);
  if( conditions.size()!=getNumberOfArguments() ) error("mismatch between number of conditions and number of arguments");

  for(unsigned i=0;i<conditions.size();++i){
     if( conditions[i]=="GT") cond[i]=true;
     else if( conditions[i]=="LT") cond[i]=false;
     else error("condition should be either GT or LT");
  }
  log.printf("  stopping calculation when value of %s is ",getPntrToArgument(0)->getName().c_str() );
  if( cond[0]==true ) log.printf("greater than ");
  else log.printf("less than ");
  log.printf("%f ",targets[0]);
  for(unsigned i=1;i<getNumberOfArguments();++i){
     log.printf("and when value of %s is ",getPntrToArgument(0)->getName().c_str() );
     if( cond[i]==true ) log.printf("greater than ");
     else log.printf("less than ");
     log.printf("%f ",targets[i]);
  }
  log.printf("\n");
  checkRead();
}

void StopWhen::update(){
  bool stopnow=true;
  for(unsigned i=0;i<getNumberOfArguments();++i){
      if( cond[i] && getArgument(i)<targets[i] ) stopnow=false;
      else if( !cond[i] && getArgument(i)>targets[i] ) stopnow=false;
  }
  if( stopnow ){
     log.printf("plumed stopped calculation as stopwhen condition %s was satisfied",getLabel().c_str());
     plumed.stop();
  }
}

}
}
