#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC FUNCTION PIECEWISE
/**

\par Examples
\verbatim
PIECEWISE  LABEL=distance  ARG=distance2 POWERS=0.5
PRINT ARG=distance,distance2
\endverbatim
(See also \ref PRINT and \ref DISTANCE).


*/
//+ENDPLUMEDOC


class FunctionPiecewise :
  public Function
{
  bool normalize;
  std::vector<std::pair<double,double> > points;
public:
  FunctionPiecewise(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(FunctionPiecewise,"PIECEWISE")

void FunctionPiecewise::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.use("PERIODIC");
  keys.add("numbered","POINT","This keyword appears multiple times as POINTx with x=0,1,2,...,n.");
  keys.reset_style("POINT","compulsory");
}

FunctionPiecewise::FunctionPiecewise(const ActionOptions&ao):
Action(ao),
Function(ao)
{
  for(int i=0;;i++){
    std::vector<double> pp;
    if(!parseNumberedVector("POINT",i,pp) ) break;
     if(pp.size()!=2) error("points should be in x,y format");
     points.push_back(std::pair<double,double>(pp[0],pp[1]));
     if(i>0 && points[i].first<=points[i-1].first) error("points abscissas should be monotonously increasing");
  }

  if(getNumberOfArguments()==1){
    addValueWithDerivatives(); 
    setNotPeriodic();
  }else{
    for(int i=0;i<getNumberOfArguments();i++){
      string s; Tools::convert(i+1,s);
      addComponentWithDerivatives(s); 
      getPntrToComponent(i)->setNotPeriodic();
    }
  }
  checkRead();

  log.printf("  on points:");
  for(unsigned i=0;i<points.size();i++) log.printf("   (%f,%f)",points[i].first,points[i].second);
  log.printf("\n");
}

void FunctionPiecewise::calculate(){
  for(int i=0;i<getNumberOfArguments();i++){
    double val=getArgument(i);
    int p=0;
    for(;p<points.size();p++){
      if(val<points[p].first) break;
    }
    double f,d;
    if(p==0){
      f=points[0].second;
      d=0.0;
    } if(p==points.size()){
      f=points[points.size()-1].second;
      d=0.0;
    } else {
      double m=(points[p].second-points[p-1].second) / (points[p].first-points[p-1].first);
      f=m*(val-points[p-1].first)+points[p-1].second;
      d=m;
    }
    if(getNumberOfArguments()==1) {
      setValue(f);
      setDerivative(i,d);
    } else {
      Value* v=getPntrToComponent(i);
      v->set(f);
      v->addDerivative(i,d);
    }
  }
}

}


