#include "Bias.h"
#include "Colvar.h"
#include <cassert>

using namespace PLMD;
using namespace std;

Bias::Bias(const ActionOptions&ao) : 
ActionWithArguments(ao),
outputForces(getNumberOfArguments(),0.0)    /// Actually this will break when readin is changed
{
}

void Bias::readBias(){
   readAction();
   std::vector<double> domain(2,0.0);
   readActionWithArguments( domain );
   addValue("energy", true, true);
   addValue("force2", true, false);
}

void Bias::apply(){
  if( onStep() ) applyForces( outputForces );
//  if(onStep()) for(unsigned i=0;i<getNumberOfArguments();++i){
//    getArguments()[i]->addForce(getStride()*outputForces[i]);
//  }
}


