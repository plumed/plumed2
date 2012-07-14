#include <cmath>
// this is an hack for cmath disambiguation

double (*cmathLog)(double) = log; 

#include "Function.h"
#include "ActionRegister.h"

#include <string>
#include <cstring>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC FUNCTION PATHMSD
/*
This file provides a template for if you want to introduce a new CV.

<!-----You should add a description of your CV here---->

\par Examples

<!---You should put an example of how to use your CV here--->

*/
//+ENDPLUMEDOC
   
class FunctionPathMSD : public Function {
  double lambda;
  bool pbc;

public:
  FunctionPathMSD(const ActionOptions&);
// active methods:
  virtual void calculate();
  static void registerKeywords(Keywords& keys);
};

PLUMED_REGISTER_ACTION(FunctionPathMSD,"PATHMSD")

void FunctionPathMSD::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
//  keys.addFlag("TEMPLATE_DEFAULT_OFF_FLAG",false,"flags that are by default not performed should be specified like this");
//  keys.addFlag("TEMPLATE_DEFAULT_ON_FLAG",true,"flags that are by default performed should be specified like this");
  keys.use("ARG");
  keys.add("compulsory","LAMBDA","all compulsory keywords should be added like this with a description here");
  keys.add("optional","NEIGHLIST","all optional keywords that have input should be added like a description here");
}

FunctionPathMSD::FunctionPathMSD(const ActionOptions&ao):
Action(ao),
Function(ao),
pbc(true)
{
//  bool nopbc=!pbc;
//  parseFlag("NOPBC",nopbc);
//  pbc=!nopbc;

  parse("LAMBDA",lambda);
  checkRead();
  log.printf("  lambda is %f\n",lambda);
  // list the action involved and check the type 
  for(unsigned i=0;i<getNumberOfArguments();i++){
       // for each value get the name and the label of the corresponding action
       std::string mylabel=getPntrToArgument(i)->getPntrToAction()->getLabel(); 
       std::string myname=getPntrToArgument(i)->getPntrToAction()->getName(); 
       log.printf("ARG  %s TYPE %s \n",mylabel.c_str(),myname.c_str());
       if(myname!="RMSD")plumed_merror("This argument is not of RMSD type!!!");
  }   
  addComponentWithDerivatives("s"); componentIsNotPeriodic("s");
  addComponentWithDerivatives("z"); componentIsNotPeriodic("z");
}
// calculator
void FunctionPathMSD::calculate(){
  double s_path=0.;
  double partition=0.;
  double tmp;
  Value* val_s_path=getPntrToComponent("s");
  Value* val_z_path=getPntrToComponent("z");
  vector<double> zvec(getNumberOfArguments());
  for(unsigned i=0;i<getNumberOfArguments();i++){ 
    zvec[i]=exp(-lambda*getArgument(i));
    s_path+=(i+1)*zvec[i];
    partition+=zvec[i];
  }
  s_path/=partition;
  val_s_path->set(s_path);
  val_z_path->set(-(1./lambda)*cmathLog(partition));
  for(unsigned i=0;i<getNumberOfArguments();i++){ 
    tmp=lambda*zvec[i]*(s_path-(i+1))/partition;
    setDerivative(val_s_path,i,tmp);
    setDerivative(val_z_path,i,zvec[i]/partition);
  }
}

}



