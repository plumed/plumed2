#ifndef __PLUMED_Bias_h
#define __PLUMED_Bias_h

#include "ActionPilot.h"
#include "ActionWithValue.h"
#include "ActionWithArguments.h"

#define PLUMED_BIAS_INIT(ao) Action(ao),Bias(ao)

namespace PLMD{

//+DEVELDOC INHERIT Bias
/**
Inherit from here if you are implementing a new form of simulation bias
*/
//+ENDDEVELDOC

class Keywords;

/// Action defining a bias which can act on other Action's
class Bias :
  public ActionPilot,
  public ActionWithValue,
  public ActionWithArguments
{
  std::vector<double> outputForces;
protected:
  void resetOutputForces();
  void setOutputForce(int i,double g);
public:
  static void registerKeywords(Keywords&);
  Bias(const ActionOptions&ao);
  void apply();
};

inline
void Bias::setOutputForce(int i,double f){
  outputForces[i]=f;
}

inline
void Bias::resetOutputForces(){
  for(unsigned i=0;i<outputForces.size();++i) outputForces[i]=0.0;
}

}

#endif

