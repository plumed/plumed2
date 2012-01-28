#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS RESTRAINT
/**
Adds harmonic and/or linear restraints on one or more variables.  

Either or both
of SLOPE and KAPPA must be present to specify the linear and harmonic force constants
respectively.  The resulting potential is given by: 
\f[
  \sum_i \frac{k_i}{2} (x_i-a_i)^2 + m_i*(x_i-a_i)
\f].

The number of components for any vector of force constants must be equal to the number
of arguments to the action.

\par Examples
The following input tells plumed to restrain the distance between atoms 3 and 5
and the distance between atoms 2 and 4, at different equilibrium
values, and to print the energy of the restraint
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
RESTRAINT ARG=d1,d2 AT=1.0,1.5 KAPPA=150.0,150.0 LABEL=restraint
PRINT ARG=restraint.bias
\endverbatim
(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC

class BiasRestraint : public Bias{
  std::vector<double> at;
  std::vector<double> kappa;
  std::vector<double> slope;
public:
  BiasRestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(BiasRestraint,"RESTRAINT")

void BiasRestraint::registerKeywords(Keywords& keys){
   Bias::registerKeywords(keys);
   keys.add("compulsory","SLOPE","0.0","specifies that the restraint is linear and what the values of the force constants on each of the variables are");
   keys.add("compulsory","KAPPA","0.0","specifies that the restraint is harmonic and what the values of the force constants on each of the variables are");
   keys.add("compulsory","AT","the position of the restraint");
}

BiasRestraint::BiasRestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
at(0),
kappa(getNumberOfArguments(),0.0),
slope(getNumberOfArguments(),0.0)
{
  parseVector("SLOPE",slope);
  assert(slope.size()==getNumberOfArguments());
  parseVector("KAPPA",kappa);
  assert(kappa.size()==getNumberOfArguments());
  parseVector("AT",at);
  assert(at.size()==getNumberOfArguments());
  checkRead();

  log.printf("  at");
  for(unsigned i=0;i<at.size();i++) log.printf(" %f",at[i]);
  log.printf("\n");
  log.printf("  with harmonic force constant");
  for(unsigned i=0;i<kappa.size();i++) log.printf(" %f",kappa[i]);
  log.printf("\n");
  log.printf("  and linear force constant");
  for(unsigned i=0;i<slope.size();i++) log.printf(" %f",slope[i]);
  log.printf("\n");

  addValue("bias");
  addValue("force2");
}


void BiasRestraint::calculate(){
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double cv=difference(i,at[i],getArgument(i));
    const double k=kappa[i];
    const double m=slope[i];
    const double f=-(k*cv+m);
    ene+=0.5*k*cv*cv+m*cv;
    setOutputForces(i,f);
    totf2+=f*f;
  };
  Value* value;
  value=getValue("bias"); setValue(value,ene);
  value=getValue("force2");  setValue(value,totf2);
}

}


