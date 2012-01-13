#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS ABMD 
/**
Adds an ratchet-and-pawl like restraint on one or more variables

\par Syntax
\verbatim
ABMD ARG=x1,x2,... KAPPA=k1,k2,... MIN=a1,a2,... TO=a1,a2,...
\endverbatim
KAPPA specifies an array of force constants, one for each variable,
and TO the target values of the restraints. Thus, the resulting potential is
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

\par Example
The following input is restraining the distance between atoms 3 and 5
and the distance between atoms 2 and 4, at different equilibrium
values, and it is printing the energy of the restraint
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
};

PLUMED_REGISTER_ACTION(BiasRatchet,"ABMD")

BiasRatchet::BiasRatchet(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
to(0),
min(getNumberOfArguments(),-1.0),
kappa(getNumberOfArguments(),0.0)
{
  parseVector("KAPPA",kappa);
  plumed_assert(kappa.size()==getNumberOfArguments());
  parseVector("MIN",min);
  plumed_assert(min.size()==getNumberOfArguments());
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

  for(unsigned i=0;i<getNumberOfArguments();i++) {char str_min[6]; sprintf(str_min,"min_%u",i+1); addValue(str_min);}
  addValue("bias");
  addValue("force2");
}


void BiasRatchet::calculate(){
  Value* value;
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0;i<getNumberOfArguments();++i){
    const double cv=difference(i,to[i],getArgument(i));
    const double cv2=cv*cv;
    const double k=kappa[i];

    if(min[i]<0.||cv2<min[i]) min[i] = cv2; 
    else {
      const double f = -2.*k*(cv2-min[i])*cv;
      setOutputForces(i,f);
      ene += 0.5*k*(cv2-min[i])*(cv2-min[i]);
      totf2+=f*f;
    }
    char str_min[6]; 
    sprintf(str_min,"min_%u",i+1); 
    value=getValue(str_min);
    setValue(value,min[i]);
  };
  value=getValue("bias"); setValue(value,ene);
  value=getValue("force2");  setValue(value,totf2);
}

}


