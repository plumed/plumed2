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
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "ReweightBase.h"

//+PLUMEDOC REWEIGHTING REWEIGHT_TEMP
/*

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace bias {

class ReweightTemperature : public ReweightBase {
private:
///
  double rtemp;
/// The biases we are using in reweighting and the args we store them separately
  std::vector<Value*> biases;
public:
  static void registerKeywords(Keywords&);
  explicit ReweightTemperature(const ActionOptions&ao);
  void prepare();
  double getLogWeight() const ;
};

PLUMED_REGISTER_ACTION(ReweightTemperature,"REWEIGHT_TEMP")

void ReweightTemperature::registerKeywords(Keywords& keys ){
  ReweightBase::registerKeywords( keys );
  keys.add("compulsory","REWEIGHT_TEMP","reweight data from a trajectory at one temperature and output the probability "
                                        "distribution at a second temperature. This is not possible during postprocessing.");
}

ReweightTemperature::ReweightTemperature(const ActionOptions&ao):
Action(ao),
ReweightBase(ao)
{
   parse("REWEIGHT_TEMP",rtemp);
   log.printf("  reweighting simulation to probabilities at temperature %f\n",rtemp);
   rtemp*=plumed.getAtoms().getKBoltzmann();

   retrieveAllBiases( "bias", biases );
}

void ReweightTemperature::prepare(){
   plumed.getAtoms().setCollectEnergy(true);
}

double ReweightTemperature::getLogWeight() const {
   // Retrieve the bias
   double bias=0.0; for(unsigned i=0;i<biases.size();++i) bias+=biases[i]->get(); 
   double energy=plumed.getAtoms().getEnergy()+bias; 
   return -( (1.0/rtemp) - (1.0/simtemp) )*(energy+bias); 
}

}
}
