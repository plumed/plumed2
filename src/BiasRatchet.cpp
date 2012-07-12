#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS ABMD 
/*
Adds a ratchet-and-pawl like restraint on one or more variables.

This action can be used to evolve a system towards a target value in
CV space using an harmonic potential moving with the thermal fluctuations of the CV
\cite ballone \cite provasi10abmd \cite camilloni11abmd. The biasing potential in this 
method is as follows:

\f$
V(\rho(t)) = \left \{ \begin{array}{ll} \frac{\alpha}{2}\left(\rho(t)-\rho_m(t)\right)^2, &\rho(t)>\rho_m(t)\\
              0, & \rho(t)\le\rho_m(t), \end{array} \right .
\f$
where
\f$
\rho(t)=\left(CV(t)-TO\right)^2
\f$
and
\f$
\rho_m(t)=\min_{0\le\tau\le t}\rho(\tau).
\f$.

The method is based on the introduction of a biasing potential which is zero when
the system is moving towards the desired arrival point and which damps the
fluctuations when the system attempts to move in the opposite direction. As in the
case of the ratchet and pawl system, propelled by thermal motion of the solvent
molecules, the biasing potential does not exert work on the system.

\par Examples
The following input sets up a restraint on the distance between atoms 3 and 5
and the distance between atoms 2 and 4, at different equilibrium
values, and tells plumed to print the energy of the restraint
\verbatim
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
ABMD ARG=d1,d2 TO=1.0,1.5 KAPPA=5.0,5.0 LABEL=abmd
PRINT ARG=abmd.bias
\endverbatim
(See also \ref DISTANCE and \ref PRINT).

*/
//+ENDPLUMEDOC

class BiasRatchet : public Bias{
  std::vector<double> to;
  std::vector<double> min;
  std::vector<double> kappa;
public:
  BiasRatchet(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(BiasRatchet,"ABMD")

void BiasRatchet::registerKeywords(Keywords& keys){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","TO","The array of target values");
  keys.add("compulsory","KAPPA","The array of force constants.");
  keys.add("optional","MIN","Array of starting values (usefull for restarting)");
}

BiasRatchet::BiasRatchet(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
to(getNumberOfArguments(),0),
min(getNumberOfArguments(),-1.0),
kappa(getNumberOfArguments(),0.0)
{
  parseVector("KAPPA",kappa);
  plumed_assert(kappa.size()==getNumberOfArguments());
  parseVector("MIN",min);
  parseVector("TO",to);
  plumed_assert(to.size()==getNumberOfArguments());
  checkRead();

  log.printf("  min");
  for(unsigned i=0;i<min.size();i++) log.printf(" %f",min[i]);
  log.printf("\n");
  log.printf("  to");
  for(unsigned i=0;i<to.size();i++) log.printf(" %f",to[i]);
  log.printf("\n");
  log.printf("  with force constant");
  for(unsigned i=0;i<kappa.size();i++) log.printf(" %f",kappa[i]);
  log.printf("\n");

  for(unsigned i=0;i<getNumberOfArguments();i++) {char str_min[6]; sprintf(str_min,"min_%u",i+1); addComponent(str_min);}
  addComponent("bias");
  addComponent("force2");
}


void BiasRatchet::calculate(){
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double cv=difference(i,to[i],getArgument(i));
    const double cv2=cv*cv;
    const double k=kappa[i];

    if(min[i]<0.||cv2<min[i]) min[i] = cv2; 
    else {
      const double f = -2.*k*(cv2-min[i])*cv;
      setOutputForce(i,f);
      ene += 0.5*k*(cv2-min[i])*(cv2-min[i]);
      totf2+=f*f;
    }
    char str_min[6]; 
    sprintf(str_min,"min_%u",i+1); 
    getPntrToComponent(str_min)->set(min[i]);
  };
  getPntrToComponent("bias")->set(ene);
  getPntrToComponent("force2")->set(totf2);
}

}


