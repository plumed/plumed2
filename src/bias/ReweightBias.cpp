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

//+PLUMEDOC REWEIGHTING REWEIGHT_BIAS
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightBias : public ReweightBase {
private:
/// The biases we are using in reweighting and the args we store them separately
  std::vector<Value*> biases;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightBias(const ActionOptions&ao);
  double getLogWeight() const ;
};

PLUMED_REGISTER_ACTION(ReweightBias,"REWEIGHT_BIAS")

void ReweightBias::registerKeywords(Keywords& keys ){
  ReweightBase::registerKeywords( keys );
}

ReweightBias::ReweightBias(const ActionOptions&ao):
Action(ao),
ReweightBase(ao)
{
   retrieveAllBiases( "bias", biases );
   if( biases.empty() ) error("there does not appear to be a bias acting on your system");
}

double ReweightBias::getLogWeight() const {
   // Retrieve the bias
   double bias=0.0; for(unsigned i=0;i<biases.size();++i) bias+=biases[i]->get();  
   return bias / simtemp;
}

}
}
