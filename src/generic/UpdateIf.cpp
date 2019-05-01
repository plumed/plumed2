/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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
#include "core/ActionPilot.h"
#include "core/ActionWithArguments.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"

using namespace std;

namespace PLMD {
namespace generic {

//+PLUMEDOC PRINTANALYSIS UPDATE_IF
/*
Conditional update of other actions.


This action can be used to enable and disable the update step for the following actions
depending on the value of its arguments. This allows for example to extract snapshots
with value of some CVs in a given range.

When called with MORE_THAN and/or LESS_THAN keywords, this action starts an if block.
The block is executed if all the arguments are less than all the respective values
in the LESS_THAN keyword (if present) and all the arguments are more than all the
respective values
in the MORE_THAN keyword (if present).

When called with the END flag, this action ends the corresponding IF block.
Notice that in this case one should also provide the ARG keyword. It is recommended to
use the same ARG keyword that was used to begin the block, so as to make the input more readable.

Of course, blocks can be nested at will.

There are many potential usages for this keyword. One might e.g. decide to analyze some variable
only when another variable is within a given range.

\warning
Notice that not all the possible usage make
particular sense. For example, conditionally updating a \ref METAD keyword
(that is: adding hills only if a variable is within a given range)
can lead to unexpected results.

\par Examples

The following input instructs plumed dump all the snapshots where an atom is in touch with
the solute.
\plumedfile
solute: GROUP ATOMS=1-124
coord: COORDINATION GROUPA=solute GROUPB=500 R_0=0.5

# A coordination number higher than 0.5 indicate that there is at least one
# atom of group `solute` at less than 5 A from atom number 500

UPDATE_IF ARG=coord MORE_THAN=0.5
DUMPATOMS ATOMS=solute,500 FILE=output.xyz
UPDATE_IF ARG=coord END
\endplumedfile

*/
//+ENDPLUMEDOC

class UpdateIf:
  public ActionPilot,
  public ActionWithArguments
{
  std::vector<double> lower;
  std::vector<double> upper;
  bool on;
  bool end;
public:
  void prepare();
  void calculate();
  void beforeUpdate();
  explicit UpdateIf(const ActionOptions&);
  static void registerKeywords(Keywords& keys);
  void apply() {}
  ~UpdateIf();
};

PLUMED_REGISTER_ACTION(UpdateIf,"UPDATE_IF")

void UpdateIf::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionPilot::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","STRIDE","1","the frequency with which the quantities of interest should be output");
  keys.addFlag("END",false,"end");
  keys.add("optional","LESS_THAN","upper bound");
  keys.add("optional","MORE_THAN","lower bound");
}

UpdateIf::UpdateIf(const ActionOptions&ao):
  Action(ao),
  ActionPilot(ao),
  ActionWithArguments(ao),
  on(false),
  end(false)
{
  parseFlag("END",end);
  parseVector("LESS_THAN",upper);
  parseVector("MORE_THAN",lower);
  if(end && upper.size()!=0) error("END and LESS_THAN are not compatible");
  if(end && lower.size()!=0) error("END and MORE_THAN are not compatible");
  if(upper.size()==0) upper.assign(getNumberOfArguments(),+std::numeric_limits<double>::max());
  if(lower.size()==0) lower.assign(getNumberOfArguments(),-std::numeric_limits<double>::max());
  if(upper.size()!=getNumberOfArguments()) error("LESS_THAN should have the same size as ARG");
  if(lower.size()!=getNumberOfArguments()) error("MORE_THAN should have the same size as ARG");
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    log<<"  boundaries for argument "<<i<<"    "<<lower[i]<<" "<<upper[i]<<"\n";
  }
  checkRead();
}

void UpdateIf::prepare() {
  on=false;
}

void UpdateIf::calculate() {
  on=true;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    if(getArgument(i)>=upper[i] || getArgument(i)<=lower[i]) on=false;
  }
}

void UpdateIf::beforeUpdate() {
  if(end) plumed.updateFlagsPop();
  else {
    if(on) plumed.updateFlagsPush(plumed.updateFlagsTop());
    else   plumed.updateFlagsPush(false);
  }
}


UpdateIf::~UpdateIf() {
}

}


}
