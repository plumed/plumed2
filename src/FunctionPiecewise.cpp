#include "ActionRegister.h"
#include "Function.h"

#include <cmath>
#include <cassert>

using namespace std;

namespace PLMD{

//+PLUMEDOC FUNCTION PIECEWISE
/**
Compute a piecewise straight line through its arguments that passes through
a set of ordered control points. For variables less than the first
(greater than the last) point, the value of the first (last) point is used.

\f[
\frac{y_{i+1}-y_i}{x_{i+1}-x_i}(s-x_i)+y_i ;  if x_i<s<x_{i+1}
\f]
\f[
y_N ; if x>x_N 
\f]
\f[
y_1 ; if x<x_1 
\f]

Control points are passed using the POINT1=... POINT2=... syntax as in the example below

If one argument is supplied, it results in a scalar quantity.
If multiple arguments are supplied, it results
in a vector of arguments.

\par Examples
\verbatim
dist1: DISTANCE ATOMS=1,10
dist2: DISTANCE ATOMS=2,11

pw: PIECEWISE POINT1=1,10 POINT2=2,pi POINT3=3,10 ARG=dist1,
ppww: PIECEWISE POINT1=1,10 POINT2=2,pi POINT3=3,10 ARG=dist1,dist2
PRINT ARG=pw,ppww.1,ppww.2
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
  keys.add("numbered","POINT","This keyword is used to specify the various points in the function above.");
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

  for(int i=0;i<getNumberOfArguments();i++)
    if(getPntrToArgument(i)->isPeriodic())
    error("Cannot use PIECEWISE on periodic arguments");

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
    } else if(p==points.size()){
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


