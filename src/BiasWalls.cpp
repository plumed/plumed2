#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS WALLS 
/**
Adds a lower/upper wall on one or more variables

\par Syntax
\verbatim
UPPER_WALLS ARG=x1,x2,... KAPPA=k1,k2,... AT=a1,a2,... OFFSET=o1,o2,... EXP=e1,e2,... EPS=s1,s2,...
LOWER_WALLS ARG=x1,x2,... KAPPA=k1,k2,... AT=a1,a2,... OFFSET=o1,o2,... EXP=e1,e2,... EPS=s1,s2,...
\endverbatim
KAPPA specifies an array of force constants, one for each variable,
and AT the center of the restraints, but with an OFFSET. Thus, the resulting potential is
\f$
  \sum_i {k_i}((x_i-a_i+o_i)/s_i)^e_i
\f$

\par Example
The following input is restraining the distance between atoms 3 and 5
and the distance between atoms 2 and 4, at different equilibrium
values, and it is printing the energy of the restraint
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
UPPER_WALLS ARG=d1,d2 AT=1.0,1.5 KAPPA=150.0,150.0 EXP=2,2 EPS=1,1 OFFSET 0,0 LABEL=uwall
LOWER_WALLS ARG=d1,d2 AT=0.0,1.0 KAPPA=150.0,150.0 EXP=2,2 EPS=1,1 OFFSET 0,0 LABEL=lwall
PRINT ARG=uwall.bias,lwall.bias
\endverbatim
(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC

class BiasUWalls : public Bias{
  std::vector<double> at;
  std::vector<double> kappa;
  std::vector<double> exp;
  std::vector<double> eps;
  std::vector<double> offset;
public:
  BiasUWalls(const ActionOptions&);
  void calculate();
};
class BiasLWalls : public Bias{
  std::vector<double> at;
  std::vector<double> kappa;
  std::vector<double> exp;
  std::vector<double> eps;
  std::vector<double> offset;
public:
  BiasLWalls(const ActionOptions&);
  void calculate();
};

PLUMED_REGISTER_ACTION(BiasUWalls,"UPPER_WALLS")
PLUMED_REGISTER_ACTION(BiasLWalls,"LOWER_WALLS")

BiasUWalls::BiasUWalls(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
at(0),
kappa(getNumberOfArguments(),0.0),
exp(getNumberOfArguments(),2.0),
eps(getNumberOfArguments(),1.0),
offset(getNumberOfArguments(),0.0)
{
  parseVector("OFFSET",kappa);
  assert(offset.size()==getNumberOfArguments());
  parseVector("EPS",kappa);
  assert(eps.size()==getNumberOfArguments());
  parseVector("EXP",kappa);
  assert(exp.size()==getNumberOfArguments());
  parseVector("KAPPA",kappa);
  assert(kappa.size()==getNumberOfArguments());
  parseVector("AT",at);
  assert(at.size()==getNumberOfArguments());
  checkRead();

  log.printf("  at");
  for(unsigned i=0;i<at.size();i++) log.printf(" %f",at[i]);
  log.printf("\n");
  log.printf("  with an offset");
  for(unsigned i=0;i<offset.size();i++) log.printf(" %f",offset[i]);
  log.printf("\n");
  log.printf("  with force constant");
  for(unsigned i=0;i<kappa.size();i++) log.printf(" %f",kappa[i]);
  log.printf("\n");
  log.printf("  and exponent");
  for(unsigned i=0;i<exp.size();i++) log.printf(" %f",exp[i]);
  log.printf("\n");
  log.printf("  rescaled");
  for(unsigned i=0;i<eps.size();i++) log.printf(" %f",eps[i]);
  log.printf("\n");

  addValue("bias");
  addValue("force2");
}

void BiasUWalls::calculate(){
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double cv=difference(i,at[i],getArgument(i));
    const double k=kappa[i];
    const double exponent=exp[i];
    const double epsilon=eps[i];
    const double off=offset[i];
    const double uscale = (cv+off)/epsilon;
    if(uscale>0.) {
      const double f=-(k/epsilon)*exponent*pow(uscale, exponent-1);
      ene+=k*pow(uscale, exponent);
      setOutputForces(i,f);
      totf2+=f*f;
    }
  };
  Value* value;
  value=getValue("bias"); setValue(value,ene);
  value=getValue("force2");  setValue(value,totf2);
}

BiasLWalls::BiasLWalls(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
at(0),
kappa(getNumberOfArguments(),0.0),
exp(getNumberOfArguments(),2.0),
eps(getNumberOfArguments(),1.0),
offset(getNumberOfArguments(),0.0)
{
  parseVector("OFFSET",kappa);
  assert(offset.size()==getNumberOfArguments());
  parseVector("EPS",kappa);
  assert(eps.size()==getNumberOfArguments());
  parseVector("EXP",kappa);
  assert(exp.size()==getNumberOfArguments());
  parseVector("KAPPA",kappa);
  assert(kappa.size()==getNumberOfArguments());
  parseVector("AT",at);
  assert(at.size()==getNumberOfArguments());
  checkRead();

  log.printf("  at");
  for(unsigned i=0;i<at.size();i++) log.printf(" %f",at[i]);
  log.printf("\n");
  log.printf("  with an offset");
  for(unsigned i=0;i<offset.size();i++) log.printf(" %f",offset[i]);
  log.printf("\n");
  log.printf("  with force constant");
  for(unsigned i=0;i<kappa.size();i++) log.printf(" %f",kappa[i]);
  log.printf("\n");
  log.printf("  and exponent");
  for(unsigned i=0;i<exp.size();i++) log.printf(" %f",exp[i]);
  log.printf("\n");
  log.printf("  rescaled");
  for(unsigned i=0;i<eps.size();i++) log.printf(" %f",eps[i]);
  log.printf("\n");

  addValue("bias");
  addValue("force2");
}

void BiasLWalls::calculate(){
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double cv=difference(i,at[i],getArgument(i));
    const double k=kappa[i];
    const double exponent=exp[i];
    const double epsilon=eps[i];
    const double off=offset[i];
    const double lscale = (cv-off)/epsilon;
    if(lscale<0.) {
      const double f=-(k/epsilon)*exponent*pow(lscale, exponent-1);
      ene+=k*pow(lscale, exponent);
      setOutputForces(i,f);
      totf2+=f*f;
    }
  };
  Value* value;
  value=getValue("bias"); setValue(value,ene);
  value=getValue("force2");  setValue(value,totf2);
}

}
