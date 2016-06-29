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
#include "core/ActionRegister.h"
#include "ReweightBase.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_METAD
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightMetad : public ReweightBase {
public:
  static void registerKeywords(Keywords&);
  explicit ReweightMetad(const ActionOptions&ao);
  double getLogWeight() const ;
};

PLUMED_REGISTER_ACTION(ReweightMetad,"REWEIGHT_METAD")

void ReweightMetad::registerKeywords(Keywords& keys ){
  ReweightBase::registerKeywords( keys ); keys.remove("ARG");
  keys.add("compulsory","ARG","*.rbias","the biases that must be taken into account when reweighting"); 
}

ReweightMetad::ReweightMetad(const ActionOptions&ao):
Action(ao),
ReweightBase(ao)
{
}

double ReweightMetad::getLogWeight() const {
   // Retrieve the bias
   double bias=0.0; for(unsigned i=0;i<getNumberOfArguments();++i) bias+=getArgument(i);  
   return bias / simtemp;
}

}
}
