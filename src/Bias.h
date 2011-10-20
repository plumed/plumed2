#ifndef __PLUMED_Bias_h
#define __PLUMED_Bias_h

#include "ActionWithArguments.h"

namespace PLMD{

/// Action defining a bias which can act on other Action's
class Bias : public ActionWithArguments {
  std::vector<double> outputForces;
protected:
  void readBias();
  void resetOutputForces();
  void setOutputForces(int i,double g);
public:
  Bias(const ActionOptions&ao);
  int setDefaultStride() const ; 
  void apply();
};

inline
int Bias::setDefaultStride() const {
  return 1;
}

inline
void Bias::setOutputForces(int i,double f){
  outputForces[i]=f;
}

inline
void Bias::resetOutputForces(){
  for(unsigned i=0;i<outputForces.size();++i) outputForces[i]=0.0;
}

}

#endif

