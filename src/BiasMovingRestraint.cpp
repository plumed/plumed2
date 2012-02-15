#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS MOVINGRESTRAINT
/**
Moving restraint (similar to \ref RESTRAINT, but dynamic)

\par Syntax
\verbatim
RESTRAINT ...
  ARG=x1,x2,...
  VERSE=v1,v2,...
  STEP0=s0 [ STEP1=s1 [STEP2=s2 [...] ] ] 
  KAPPA0=k01,k02,... [ KAPPA1=k11,k12,... [ KAPPA2=k21,k22,... [...] ] ] 
  AT0=a01,a02,...    [ AT1=a11,a12,...    [ AT1=a21,a22,...    [...] ] ]
... RESTRAINT
\endverbatim
The STEPx keyword, with x=0,1,2,...,n respresent the MD step at
which the restraint parameters take value KAPPAx and ATx. The meaning
of KAPPA and AT is the same as for \ref RESTRAINT. For step numbers less than
STEP0 or larger than STEPn, parameters for x=0 and x=n are used respectively.
For intermediate steps, parameters are linearly interpolated. If
a parameter is missing, its value at the present step is kept.
The VERSE keyword can set for each variable is the restraint is
only acting for CV larger ("U") or smaller ("L") than the
restraint. Default is "B" which is acting in both cases.

\par Examples
The following input is dragging the distance between atoms 2 and 4
from 1 to 2 in the first 1000 steps, then back in the next 1000 steps.
In the following 500 steps the restraint is progressively switched off.
\verbatim
DISTANCE ATOMS=2,4 LABEL=d
MOVINGRESTRAINT ...
  ARG=d1,d2
  STEP0=0    AT0=1.0 KAPPA0=100.0
  STEP1=1000 AT1=2.0
  STEP2=2000 AT2=1.0
  STEP3=2500         KAPPA3=0.0
... MOVINGRESTRAINT
\endverbatim
The following input is progressively building restraints
distances between atoms 1 and 5 and between atoms 2 and 4
in the first 1000 steps. Afterwards, the restraint is kept
static.
\verbatim
DISTANCE ATOMS=1,5 LABEL=d1
DISTANCE ATOMS=2,4 LABEL=d2
MOVINGRESTRAINT ...
  ARG=d1,d2 
  STEP0=0    AT0=1.0,1.5 KAPPA0=0.0,0.0
  STEP0=1000 AT0=1.0,1.5 KAPPA0=1.0,1.0
... MOVINGRESTRAINT
\endverbatim
The following input is progressively bringing atoms 1 and 2
close to each other with an upper wall
\verbatim
DISTANCE ATOMS=1,2 LABEL=d1
MOVINGRESTRAINT ...
  ARG=d1
  STEP0=0    AT0=1.0 KAPPA0=10.0
  STEP0=1000 AT0=0.0
... MOVINGRESTRAINT
\endverbatim


(See also \ref DISTANCE).

*/
//+ENDPLUMEDOC


class BiasMovingRestraint : public Bias{
  std::vector<std::vector<double> > at;
  std::vector<std::vector<double> > kappa;
  std::vector<int> step;
  vector<string> verse;
public:
  BiasMovingRestraint(const ActionOptions&);
  void calculate();
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(BiasMovingRestraint,"MOVINGRESTRAINT")

void BiasMovingRestraint::registerKeywords( Keywords& keys ){
  Bias::registerKeywords(keys);
  keys.use("ARG");
  keys.add("compulsory","VERSE","B","Tells plumed whether the restraint is only acting for CV larger (U) or smaller (L) than the restraint or whether it is acting on both sides (B)");
  keys.add("numbered","STEP","This keyword appears multiple times as STEPx with x=0,1,2,...,n.  Each value given represents the MD step at which the restraint parameters take the values KAPPAx and ATx.");
  keys.add("numbered","AT","ATx is equal to the position of the restraint at time STEPx.  For intermediate times this parameter is linearly interpolated.  If no ATx is specified for STEPx then the values of AT are kept constant during the interval of time between STEPx-1 and STEPx.");
  keys.add("numbered","KAPPA","KAPPAx is equal to the value of the force constants at time STEPx.  For intermediate times this parameter is linearly interpolated.  If no ATx is specified for STEPx then the values of KAPPAx are kept constant during the interval of time between STEPx-1 and STEPx.");
}

BiasMovingRestraint::BiasMovingRestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao),
verse(getNumberOfArguments())
{
  parseVector("VERSE",verse);
  vector<int> ss(1); ss[0]=-1;
  std::vector<double> kk( getNumberOfArguments() ), aa( getNumberOfArguments() );
  for(int i=0;;i++){
    // Read in step 
    if( !parseNumberedVector("STEP",i,ss) ) break;
    for(unsigned j=0;j<step.size();j++) plumed_massert(ss[0]>step[j],"in moving restraint step number must always increase");
    step.push_back(ss[0]);

    // Try to read kappa
    if( !parseNumberedVector("KAPPA",i,kk) ) kk=kappa[i-1]; 
    kappa.push_back(kk);

    // Now read AT
    if( !parseNumberedVector("AT",i,aa) ) aa=at[i-1];
    at.push_back(aa);
  }
  checkRead();

  for(unsigned i=0;i<step.size();i++){
    log.printf("  step%d %d\n",i,step[i]);
    log.printf("  at");
    for(unsigned j=0;j<at[i].size();j++) log.printf(" %f",at[i][j]);
    log.printf("\n");
    log.printf("  with force constant");
    for(unsigned j=0;j<kappa[i].size();j++) log.printf(" %f",kappa[i][j]);
    log.printf("\n");
  };

  addComponent("bias");
  addComponent("force2");
  addComponent("work");
}


void BiasMovingRestraint::calculate(){
  double ene=0.0;
  double totf2=0.0;
  unsigned narg=getNumberOfArguments();
  int now=getStep();
  std::vector<double> kk(narg),aa(narg);
  if(now<=step[0]){
    kk=kappa[0];
    aa=at[0];
  } else if(now>=step[step.size()-1]){
    kk=kappa[step.size()-1];
    aa=at[step.size()-1];
  } else {
    unsigned i=0;
    for(i=1;i<step.size();i++) if(now<step[i]) break;
    double c2=(now-step[i-1])/double(step[i]-step[i-1]);
    double c1=1.0-c2;
    for(unsigned j=0;j<narg;j++) kk[j]=(c1*kappa[i-1][j]+c2*kappa[i][j]);
    for(unsigned j=0;j<narg;j++) aa[j]=(c1*at[i-1][j]+c2*at[i][j]);
  }
  for(unsigned i=0;i<narg;++i){
    const double cv=difference(i,aa[i],getArgument(i));
    const double k=kk[i];
    const double f=-k*cv;
    if(verse[i]=="U" && cv<0) continue;
    if(verse[i]=="L" && cv>0) continue;
    assert(verse[i]=="U" || verse[i]=="L" || verse[i]=="B");
    ene+=0.5*k*cv*cv;
    setOutputForce(i,f);
    totf2+=f*f;
  };
  getPntrToComponent("bias")->set(ene);
  getPntrToComponent("force2")->set(totf2);
}

}


