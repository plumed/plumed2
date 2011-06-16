#include "Bias.h"
#include "ActionRegister.h"

#include <cassert>

using namespace std;


namespace PLMD{

//+PLUMEDOC BIAS MOVINGRESTRAINT
/**
Moving restraint

Similar to \ref RESTRAINT but can be moved during the simulation
   
\par Syntax

Example
\verbatim
MOVINGRESTRAINT
\endverbatim

If you use this variable, please cite the following work ... ... ...

*/
//+ENDPLUMEDOC


class BiasMovingRestraint : public Bias{
  std::vector<std::vector<double> > at;
  std::vector<std::vector<double> > kappa;
  std::vector<int> step;
public:
  BiasMovingRestraint(const ActionOptions&);
  void calculate();
};

PLUMED_REGISTER_ACTION(BiasMovingRestraint,"MOVINGRESTRAINT")

BiasMovingRestraint::BiasMovingRestraint(const ActionOptions&ao):
PLUMED_BIAS_INIT(ao)
{
  for(int i=0;;i++){
    int ss=-1;
    std::vector<double> kk,aa;
    string s;
    Tools::convert(i,s);
    parse("STEP"+s,ss);
    if(ss<0) break;
    for(unsigned j=0;j<step.size();j++) assert(ss>step[j]);
    step.push_back(ss);
    parseVector("KAPPA"+s,kk);
    if(kk.size()==0 && i>0) kk=kappa[i-1];
    assert(kk.size()==getNumberOfArguments());
    kappa.push_back(kk);
    parseVector("AT"+s,aa);
    if(at.size()==0 && i>0) aa=at[i-1];
    assert(aa.size()==getNumberOfArguments());
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

  addValue("");
  addValue("Energy");
  addValue("Force2");
  addValue("Work");
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
    const double cv=(getArgument(i)-aa[i]);
    const double k=kk[i];
    const double f=-k*cv;
    ene+=0.5*k*cv*cv;
    setOutputForces(i,f);
    totf2+=f*f;
  };
  getValue("")->setValue(ene);
  getValue("Energy")->setValue(ene);
  getValue("Force2")->setValue(totf2);
}

}


