/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The plumed team
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
#include "ReweightBase.h"
#include "core/ActionRegister.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_WELLTEMPERED
/*

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightWellTempered : public ReweightBase {
private:
  double logheight, deltaT;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightWellTempered(const ActionOptions&ao);
  double getLogWeight();
};

PLUMED_REGISTER_ACTION(ReweightWellTempered,"REWEIGHT_WELLTEMPERED")

void ReweightWellTempered::registerKeywords(Keywords& keys ) {
  ReweightBase::registerKeywords( keys ); keys.remove("ARG");
  keys.add("optional","ARG","the biases that must be taken into account when reweighting");
  keys.add("compulsory","BIASFACTOR","the well tempered factor");
  keys.add("compulsory","HEIGHT","the initial heights of the Guassians");
}

ReweightWellTempered::ReweightWellTempered(const ActionOptions&ao):
  Action(ao),
  ReweightBase(ao)
{
  double biasf, height; parse("BIASFACTOR",biasf); parse("HEIGHT",height);
  log.printf("  initial weights are %f and bias factor is set equal to %f \n", height, biasf );
  deltaT = simtemp*(biasf - 1.0); logheight = std::log( height );
}

double ReweightWellTempered::getLogWeight() {
  if( getNumberOfArguments()==0 ) error("bias was not set in input");
  // Retrieve the bias
  double bias=0.0; for(unsigned i=0; i<getNumberOfArguments(); ++i) bias+=getArgumentScalar(i);
  return -bias / deltaT + logheight;
}

}
}
