/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
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
#include "Bias.h"
#include "tools/Random.h"
#include "ActionRegister.h"
#include <ctime>

namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS ABMD
/*
Adds a ratchet-and-pawl like restraint on one or more variables.

This action can be used to evolve a system towards a target value in
CV space using an harmonic potential moving with the thermal fluctuations of the CV
\cite ballone \cite provasi10abmd \cite camilloni11abmd. The biasing potential in this
method is as follows:

\f$
V(\rho(t)) = \left \{ \begin{array}{ll} \frac{K}{2}\left(\rho(t)-\rho_m(t)\right)^2, &\rho(t)>\rho_m(t)\\
              0, & \rho(t)\le\rho_m(t), \end{array} \right .
\f$


where


\f$
\rho(t)=\left(CV(t)-TO\right)^2
\f$


and


\f$
\rho_m(t)=\min_{0\le\tau\le t}\rho(\tau)+\eta(t)
\f$.

The method is based on the introduction of a biasing potential which is zero when
the system is moving towards the desired arrival point and which damps the
fluctuations when the system attempts to move in the opposite direction. As in the
case of the ratchet and pawl system, propelled by thermal motion of the solvent
molecules, the biasing potential does not exert work on the system. \f$\eta(t)\f$ is
an additional white noise acting on the minimum position of the bias.

\par Examples

The following input sets up two biases, one on the distance between atoms 3 and 5
and another on the distance between atoms 2 and 4. The two target values are defined
using TO and the two strength using KAPPA. The total energy of the bias is printed.
\plumedfile
DISTANCE ATOMS=3,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
ABMD ARG=d1,d2 TO=1.0,1.5 KAPPA=5.0,5.0 LABEL=abmd
PRINT ARG=abmd.bias,abmd.d1_min,abmd.d2_min
\endplumedfile

*/
//+ENDPLUMEDOC

class ABMD : public Bias {
  std::vector<double> to;
  std::vector<double> min;
  std::vector<double> kappa;
  std::vector<double> temp;
  std::vector<int> seed;
  std::vector<Random> random;
public:
  explicit ABMD(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(ABMD,"ABMD")

void ABMD::registerKeywords(Keywords& keys) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","TO","The array of target values");
  keys.add("compulsory","KAPPA","The array of force constants.");
  keys.add("optional","MIN","Array of starting values for the bias (set rho_m(t), otherwise it is set using the current value of ARG)");
  keys.add("optional","NOISE","Array of white noise intensities (add a temperature to the ABMD)");
  keys.add("optional","SEED","Array of seeds for the white noise (add a temperature to the ABMD)");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
  keys.addOutputComponent("_min","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          " These quantities will be named with the arguments of the bias followed by "
                          "the character string _min. These quantities tell the user the minimum value assumed by rho_m(t).");
}

ABMD::ABMD(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  to(getNumberOfArguments(),0.0),
  min(getNumberOfArguments(),-1.0),
  kappa(getNumberOfArguments(),0.0),
  temp(getNumberOfArguments(),0.0),
  seed(getNumberOfArguments(),std::time(0)),
  random(getNumberOfArguments())
{
  // Note : parseVector will check that number of arguments is correct
  parseVector("KAPPA",kappa);
  parseVector("MIN",min);
  if(min.size()==0) min.assign(getNumberOfArguments(),-1.0);
  if(min.size()!=getNumberOfArguments()) error("MIN array should have the same size as ARG array");
  parseVector("NOISE",temp);
  parseVector("SEED",seed);
  parseVector("TO",to);
  checkRead();

  log.printf("  min");
  for(unsigned i=0; i<min.size(); i++) log.printf(" %f",min[i]);
  log.printf("\n");
  log.printf("  to");
  for(unsigned i=0; i<to.size(); i++) log.printf(" %f",to[i]);
  log.printf("\n");
  log.printf("  with force constant");
  for(unsigned i=0; i<kappa.size(); i++) log.printf(" %f",kappa[i]);
  log.printf("\n");

  for(unsigned i=0; i<getNumberOfArguments(); i++) {
    std::string str_min=getPntrToArgument(i)->getName()+"_min";
    addComponent(str_min); componentIsNotPeriodic(str_min);
    if(min[i]!=-1.0) getPntrToComponent(str_min)->set(min[i]);
  }
  for(unsigned i=0; i<getNumberOfArguments(); i++) {random[i].setSeed(-seed[i]);}
  addComponent("force2"); componentIsNotPeriodic("force2");
}


void ABMD::calculate() {
  double ene=0.0;
  double totf2=0.0;
  for(unsigned i=0; i<getNumberOfArguments(); ++i) {
    const double cv=difference(i,to[i],getArgument(i));
    const double cv2=cv*cv;
    const double k=kappa[i];
    double noise=0.;
    double diff=temp[i];
    if(diff>0) {
      noise = 2.*random[i].Gaussian()*diff;
      if(cv2<=diff) { diff=0.; temp[i]=0.; }
    }

    // min < 0 means that the variable has not been used in the input file, so the current position of the CV is used
    // cv2 < min means that the collective variable is nearer to the target value than at any other previous time so
    // min is set to the CV value
    if(min[i]<0.||cv2<min[i]) {
      min[i] = cv2;
    } else {
      // otherwise a noise is added to the minimum value
      min[i] += noise;
      const double f = -2.*k*(cv2-min[i])*cv;
      setOutputForce(i,f);
      ene += 0.5*k*(cv2-min[i])*(cv2-min[i]);
      totf2+=f*f;
    }
    getPntrToComponent(i+1)->set(min[i]);
  }
  setBias(ene);
  getPntrToComponent("force2")->set(totf2);
}

}
}


