#include "Value.h"
#include "ActionWithValue.h"

using namespace PLMD;

Value::Value( ActionWithValue&action,const std::string& name, const unsigned& nd, const std::vector<double>& domain ):
action(action),
value(0.0),
myname(name),
deriv(false),
periodicity(unset),
min(0.0),
max(0.0)
{
  assert(domain.size()==2); 
  if( domain[0]==0.0 && domain[1]==0.0 ){
     periodicity=notperiodic;
  } else {
     periodicity=periodic;
     min=domain[0]; max=domain[1];
  } 
  derivatives.resize(nd);
}

double Value::difference(double d1,double d2) const {
  assert(periodicity!=unset);
  if(periodicity==periodic){
    double s=(d2-d1)/(max-min);
    s=Tools::pbc(s);
    return s*(max-min);
  } else{
    return d2-d1;
  }
} 
