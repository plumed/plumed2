#include "Function.h"
#include "Colvar.h"
#include <cassert>

using namespace PLMD;
using namespace std;

Function::Function(const ActionOptions&ao) : ActionWithArguments(ao) {
  forbidKeyword("STRIDE");
  registerKeyword(1, "PERIODIC", "if the function is not periodic should be PERIODIC=NO otherwise specfies the domain"); 
}

void Function::readFunction(){
  readAction();
  std::vector<double> domain(2,0.0);
  std::vector<std::string> period; parseVector("PERIODIC",period);
  if(period.size()==1 && period[0]=="NO"){
     readActionWithArguments( domain );
  } else if(period.size()==2 && Tools::convert(period[0],domain[0]) && Tools::convert(period[1],domain[1])){
     readActionWithArguments( domain );  
  } else {
     error("input to PERIODIC keyword makes no sense should be NO or the domain of the function");
  }
}

void Function::apply(){

  vector<double> f( getNumberOfArguments(), 0.0 ), forces( getNumberOfArguments(),0.0 );
  bool at_least_one_forced=false;

  for(unsigned i=0;i<getNumberOfValues();++i){
    if( getForces( i, forces ) ){
       at_least_one_forced=true;
       for(unsigned j=0;j<forces.size();j++){ f[j]+=forces[j]; }
    } 
    //if(!getValue(i)->checkForced()) continue;
    //at_least_one_forced=true;
    //const vector<double> & derivatives(getValue(i)->getDerivatives());
    //for(unsigned j=0;j<derivatives.size();j++){
    //  f[j]+=getForce(i)*derivatives[j];
    //}
  }
  if(at_least_one_forced) applyForces( f );
//     for(unsigned i=0;i<getNumberOfArguments();++i){
//       //getArguments()[i]->addForce(f[i]);
//     }
//  }
}



