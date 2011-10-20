#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS MOVINGRESTRAINT
/**
Moving restraint (similar to \ref RESTRAINT, but dynamic)

\par Example
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
};

PLUMED_REGISTER_ACTION(BiasMovingRestraint,"MOVINGRESTRAINT")

BiasMovingRestraint::BiasMovingRestraint(const ActionOptions&ao):
Bias(ao)
{
  // GAT register keywords
  registerKeyword(3,"STEP","In a moving restraint the bias changes as a function of time.  These changes are defined using instances of the STEP keyword - i.e. STEP1, STEP2, ....  In essence this keyword states that at the specified time the value of kappa should be equal to the corresponding value of kappa and the value of at should be equal to the corresponding value of at.  That is to say at the time specified by STEP1, kappa should equal KAPPA1 and at should equal AT1.  Inbetween the times specified by STEP keywords the force constants and restraint positions are linearly interpolated");
  registerKeyword(3,"KAPPA","The values of the force constant at the times specified by STEP during the simulation");
  registerKeyword(3,"AT","The location of the center of the restraint at the times specified the STEP during the simulation");
  registerKeyword(0,"VERSE","(default=B) how is the bias enforcing the restraint.  U indicates that it only acts if the instantaneous value of the CV is larger than the current location of the restraint, L indicates that it only acts if the instantaneous value of the CV is less than the current location of the restraint and B indicates that the restraint is always acting.");
  readBias();
  
  parseVector("VERSE",verse);
  if(verse.size()==0) verse.assign(getNumberOfArguments(),"B");

  if( verse.size()!=getNumberOfArguments() ) error("wrong number of arguments in verse keyword");
  for(unsigned i=0;i<verse.size();++i){
     if( verse[i]!="U" && verse[i]!="L" && verse[i]!="B" ) error( verse[i] + " is not a valid input for the VERSE keyword.  Use U, L or B");
  }

  for(int i=0;;i++){
    //int ss=-1;
    std::vector<int> ss;
    std::vector<double> kk,aa;
    string s; Tools::convert(i,s);
    if( !parseNumberedVector("STEP", i, ss) ) break;
    if( ss.size()!=1 ) error("argument to step keyword should be a single number"); 
    //parse("STEP"+s,ss);
    //if(ss<0) break;
    for(unsigned j=0;j<step.size();j++){
      // assert(ss>step[j]);
      if( ss[0]<step[j] ) error("really - you want you're simulation to go backwards in time!");
    }
    step.push_back(ss[0]);
    if( !parseNumberedVector("KAPPA", i, kk) ) error("for STEP" + s + " there is no corresponding force constant specified use KAPPA" + s);
    if(kk.size()==0 && i>0) kk=kappa[i-1];
    if( kk.size()!=getNumberOfArguments() ) error("force contant " + s + " has wrong size");
    //assert(kk.size()==getNumberOfArguments());
    kappa.push_back(kk);
    if( !parseNumberedVector("AT", i, aa) ) error("for STEP" + s + " there is no corresponding restraint position specified use AT" + s);
    if(aa.size()==0 && i>0) aa=at[i-1];
    if ( aa.size()!=getNumberOfArguments() ) error("restraint location " + s + " has wrong size");
    //assert(aa.size()==getNumberOfArguments());
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
  addValue("Work", true, false );
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
    //assert(verse[i]=="U" || verse[i]=="L" || verse[i]=="B");
    ene+=0.5*k*cv*cv;
    setOutputForces(i,f);
    totf2+=f*f;
  };
  //Value* value;
  //value=getValue("Energy");
  //setValue(value,ene);
  //value=getValue("Force2");
  //setValue(value,totf2);
}

}


