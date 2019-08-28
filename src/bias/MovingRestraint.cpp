/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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
#include "ActionRegister.h"


using namespace std;


namespace PLMD {
namespace bias {

//+PLUMEDOC BIAS MOVINGRESTRAINT
/*
Add a time-dependent, harmonic restraint on one or more variables.

This form of bias can be used to performed steered MD \cite Grubmuller3
and Jarzynski sampling \cite jarzynski.

The harmonic restraint on your system is given by:

\f[
V(\vec{s},t) = \frac{1}{2} \kappa(t) ( \vec{s} - \vec{s}_0(t) )^2
\f]

The time dependence of \f$\kappa\f$ and \f$\vec{s}_0\f$ are specified by a list of
STEP, KAPPA and AT keywords.  These keywords tell plumed what values \f$\kappa\f$ and \f$\vec{s}_0\f$
should have at the time specified by the corresponding STEP keyword.  In between these times
the values of \f$\kappa\f$ and \f$\vec{s}_0\f$ are linearly interpolated.

Additional material and examples can be also found in the tutorial \ref belfast-5

\par Examples

The following input is dragging the distance between atoms 2 and 4
from 1 to 2 in the first 1000 steps, then back in the next 1000 steps.
In the following 500 steps the restraint is progressively switched off.
\plumedfile
DISTANCE ATOMS=2,4 LABEL=d
MOVINGRESTRAINT ...
  ARG=d
  STEP0=0    AT0=1.0 KAPPA0=100.0
  STEP1=1000 AT1=2.0
  STEP2=2000 AT2=1.0
  STEP3=2500         KAPPA3=0.0
... MOVINGRESTRAINT
\endplumedfile
The following input is progressively building restraints
distances between atoms 1 and 5 and between atoms 2 and 4
in the first 1000 steps. Afterwards, the restraint is kept
static.
\plumedfile
DISTANCE ATOMS=1,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
MOVINGRESTRAINT ...
  ARG=d1,d2
  STEP0=0    AT0=1.0,1.5 KAPPA0=0.0,0.0
  STEP1=1000 AT1=1.0,1.5 KAPPA1=1.0,1.0
... MOVINGRESTRAINT
\endplumedfile
The following input is progressively bringing atoms 1 and 2
close to each other with an upper wall
\plumedfile
DISTANCE ATOMS=1,2 LABEL=d1
MOVINGRESTRAINT ...
  ARG=d1
  VERSE=U
  STEP0=0    AT0=1.0 KAPPA0=10.0
  STEP1=1000 AT1=0.0
... MOVINGRESTRAINT
\endplumedfile

By default the Action is issuing some values which are
the work on each degree of freedom, the center of the harmonic potential,
the total bias deposited

(See also \ref DISTANCE).

\attention Work is not computed properly when KAPPA is time dependent.

*/
//+ENDPLUMEDOC


class MovingRestraint : public Bias {
  std::vector<std::vector<double> > at;
  std::vector<std::vector<double> > kappa;
  std::vector<long int> step;
  std::vector<double> oldaa;
  std::vector<double> oldk;
  std::vector<double> olddpotdk;
  std::vector<double> oldf;
  std::vector<string> verse;
  std::vector<double> work;
  double tot_work;
public:
  explicit MovingRestraint(const ActionOptions&);
  void calculate() override;
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(MovingRestraint,"MOVINGRESTRAINT")

void MovingRestraint::registerKeywords( Keywords& keys ) {
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","VERSE","B","Tells plumed whether the restraint is only acting for CV larger (U) or smaller (L) than "
           "the restraint or whether it is acting on both sides (B)");
  keys.add("numbered","STEP","This keyword appears multiple times as STEP\\f$x\\f$ with x=0,1,2,...,n. Each value given represents "
           "the MD step at which the restraint parameters take the values KAPPA\\f$x\\f$ and AT\\f$x\\f$.");
  keys.reset_style("STEP","compulsory");
  keys.add("numbered","AT","AT\\f$x\\f$ is equal to the position of the restraint at time STEP\\f$x\\f$. For intermediate times this parameter "
           "is linearly interpolated. If no AT\\f$x\\f$ is specified for STEP\\f$x\\f$ then the values of AT are kept constant "
           "during the interval of time between STEP\\f$x-1\\f$ and STEP\\f$x\\f$.");
  keys.reset_style("AT","compulsory");
  keys.add("numbered","KAPPA","KAPPA\\f$x\\f$ is equal to the value of the force constants at time STEP\\f$x\\f$. For intermediate times this "
           "parameter is linearly interpolated.  If no KAPPA\\f$x\\f$ is specified for STEP\\f$x\\f$ then the values of KAPPA\\f$x\\f$ "
           "are kept constant during the interval of time between STEP\\f$x-1\\f$ and STEP\\f$x\\f$.");
  keys.reset_style("KAPPA","compulsory");
  keys.addOutputComponent("work","default","the total work performed changing this restraint");
  keys.addOutputComponent("force2","default","the instantaneous value of the squared force due to this bias potential");
  keys.addOutputComponent("_cntr","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "these quantities will named with  the arguments of the bias followed by "
                          "the character string _cntr. These quantities give the instantaneous position "
                          "of the center of the harmonic potential.");
  keys.addOutputComponent("_work","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "These quantities will named with the arguments of the bias followed by "
                          "the character string _work. These quantities tell the user how much work has "
                          "been done by the potential in dragging the system along the various colvar axis.");
  keys.addOutputComponent("_kappa","default","one or multiple instances of this quantity can be referenced elsewhere in the input file. "
                          "These quantities will named with the arguments of the bias followed by "
                          "the character string _kappa. These quantities tell the user the time dependent value of kappa.");
}

MovingRestraint::MovingRestraint(const ActionOptions&ao):
  PLUMED_BIAS_INIT(ao),
  verse(getNumberOfArguments())
{
  parseVector("VERSE",verse);
  vector<long int> ss(1); ss[0]=-1;
  std::vector<double> kk( getNumberOfArguments() ), aa( getNumberOfArguments() );
  for(int i=0;; i++) {
    // Read in step
    if( !parseNumberedVector("STEP",i,ss) ) break;
    for(unsigned j=0; j<step.size(); j++) {
      if(ss[0]<step[j]) error("in moving restraint step number must always increase");
    }
    step.push_back(ss[0]);

    // Try to read kappa
    if( !parseNumberedVector("KAPPA",i,kk) ) kk=kappa[i-1];
    kappa.push_back(kk);

    // Now read AT
    if( !parseNumberedVector("AT",i,aa) ) aa=at[i-1];
    at.push_back(aa);
  }
  checkRead();

  for(unsigned i=0; i<step.size(); i++) {
    log.printf("  step%u %ld\n",i,step[i]);
    log.printf("  at");
    for(unsigned j=0; j<at[i].size(); j++) log.printf(" %f",at[i][j]);
    log.printf("\n");
    log.printf("  with force constant");
    for(unsigned j=0; j<kappa[i].size(); j++) log.printf(" %f",kappa[i][j]);
    log.printf("\n");
  };

  addComponent("force2"); componentIsNotPeriodic("force2");

  // add the centers of the restraint as additional components that can be retrieved (useful for debug)

  std::string comp;
  for(unsigned i=0; i< getNumberOfArguments() ; i++) {
    comp=getPntrToArgument(i)->getName()+"_cntr"; // each spring has its own center
    addComponent(comp); componentIsNotPeriodic(comp);
    comp=getPntrToArgument(i)->getName()+"_work"; // each spring has its own work
    addComponent(comp); componentIsNotPeriodic(comp);
    comp=getPntrToArgument(i)->getName()+"_kappa"; // each spring has its own kappa
    addComponent(comp); componentIsNotPeriodic(comp);
    work.push_back(0.); // initialize the work value
  }
  addComponent("work"); componentIsNotPeriodic("work");
  tot_work=0.0;

  log<<"  Bibliography ";
  log<<cite("Grubmuller, Heymann, and Tavan, Science 271, 997 (1996)")<<"\n";

}


void MovingRestraint::calculate() {
  double ene=0.0;
  double totf2=0.0;
  unsigned narg=getNumberOfArguments();
  long int now=getStep();
  std::vector<double> kk(narg),aa(narg),f(narg),dpotdk(narg);
  if(now<=step[0]) {
    kk=kappa[0];
    aa=at[0];
    oldaa=at[0];
    oldk=kappa[0];
    olddpotdk.resize(narg);
    oldf.resize(narg);
  } else if(now>=step[step.size()-1]) {
    kk=kappa[step.size()-1];
    aa=at[step.size()-1];
  } else {
    unsigned i=0;
    for(i=1; i<step.size(); i++) if(now<step[i]) break;
    double c2=(now-step[i-1])/double(step[i]-step[i-1]);
    double c1=1.0-c2;
    for(unsigned j=0; j<narg; j++) kk[j]=(c1*kappa[i-1][j]+c2*kappa[i][j]);
    for(unsigned j=0; j<narg; j++) {
      const double a1=at[i-1][j];
      const double a2=at[i][j];
      aa[j]=(c1*a1+c2*(a1+difference(j,a1,a2)));
    }
  }
  tot_work=0.0;
  for(unsigned i=0; i<narg; ++i) {
    const double cv=difference(i,aa[i],getArgument(i)); // this gives: getArgument(i) - aa[i]
    getPntrToComponent(getPntrToArgument(i)->getName()+"_cntr")->set(aa[i]);
    const double k=kk[i];
    f[i]=-k*cv;
    if(verse[i]=="U" && cv<0) continue;
    if(verse[i]=="L" && cv>0) continue;
    plumed_assert(verse[i]=="U" || verse[i]=="L" || verse[i]=="B");
    dpotdk[i]=0.5*cv*cv;
    if(oldaa.size()==aa.size() && oldf.size()==f.size()) work[i]+=0.5*(oldf[i]+f[i])*(aa[i]-oldaa[i]) + 0.5*( dpotdk[i]+olddpotdk[i] )*(kk[i]-oldk[i]);
    getPntrToComponent(getPntrToArgument(i)->getName()+"_work")->set(work[i]);
    getPntrToComponent(getPntrToArgument(i)->getName()+"_kappa")->set(kk[i]);
    tot_work+=work[i];
    ene+=0.5*k*cv*cv;
    setOutputForce(i,f[i]);
    totf2+=f[i]*f[i];
  };
  getPntrToComponent("work")->set(tot_work);
  oldf=f;
  oldaa=aa;
  oldk=kk;
  olddpotdk=dpotdk;
  setBias(ene);
  getPntrToComponent("force2")->set(totf2);
}

}
}


