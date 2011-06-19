#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS RESTRAINT
/**
Harmonic restraint

This bias imposes an harmonic restraint on one or more
collective variables (or function). Besides the
usual STRIDE, ARG and LABEL, it accepts the following
keywords:
- KAPPA an array of force constants (in PLUMED units)
- AT an array defining the center of the restraint
Example
\verbatim
DISTANCE ATOMS=5,10 LABEL=d1
DISTANCE ATOMS=6,11 LABEL=d2
DISTANCE ATOMS=7,12 LABEL=d3
RESTRAINT ARG=d1,d2,d3 KAPPA=0.1,0.1,0.5 AT=5.0,6.0,7.0
\endverbatim

*/
//+ENDPLUMEDOC

class BiasRestraint : public Bias{
  std::vector<double> at;
  std::vector<double> kappa;
public:
  BiasRestraint(const ActionOptions&);
  void calculate();
};

PLUMED_REGISTER_ACTION(BiasRestraint,"RESTRAINT")

BiasRestraint::BiasRestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
at(0),
kappa(getNumberOfArguments(),0.0)
{
  parseVector("KAPPA",kappa);
  assert(kappa.size()==getNumberOfArguments());
  parseVector("AT",at);
  assert(at.size()==getNumberOfArguments());
  checkRead();

  log.printf("  at");
  for(unsigned i=0;i<at.size();i++) log.printf(" %f",at[i]);
  log.printf("\n");
  log.printf("  with force constant");
  for(unsigned i=0;i<kappa.size();i++) log.printf(" %f",kappa[i]);
  log.printf("\n");

  addValue("");
  addValue("Energy");
  addValue("Force2");
}


void BiasRestraint::calculate(){
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double cv=difference(i,at[i],getArgument(i));
    const double k=kappa[i];
    const double f=-k*cv;
    ene+=0.5*k*cv*cv;
    setOutputForces(i,f);
    totf2+=f*f;
  };
  setValue(ene);
  Value* value;
  value=getValue("Energy"); setValue(value,ene);
  value=getValue("Force2");  setValue(value,totf2);
}

}


