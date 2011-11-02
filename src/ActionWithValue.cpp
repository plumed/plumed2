#include "ActionWithValue.h"
#include "PlumedMain.h"

using namespace std;
using namespace PLMD;

ActionWithValue::ActionWithValue(const ActionOptions&ao):
  Action(ao),
  numericalDerivatives(false),
  parallel(false)
{
  registerKeyword(0, "NUMERICAL_DERIVATIVES", "calculate the derivatives for these quantities numerically"); 
  registerKeyword(0, "PARALLELIZE", "(default=off) if many quantities have to be calculated to get your colvar then you can parallize the calculation using this keyword.  For example, if you are calculating the number of segements of the backbone that are within 1.0 A RMSD of a perfect alpha helix using ALPHARMSD LESS_THAN=0.1 BACKBONE=1-50 this calculation involves many RMSD calculations.  It may therefore be worthwhile to parallize the process using this keyword." );
}

ActionWithValue::~ActionWithValue(){
  for(unsigned i=0;i<values.size();++i) delete values[i];
}

void ActionWithValue::readActionWithValue( const unsigned& nd, const std::vector<double>& d ){
  parseFlag("NUMERICAL_DERIVATIVES",numericalDerivatives); 
  parseFlag("PARALLELIZE",parallel);
  if (numericalDerivatives) log.printf("  using numerical derivatives\n");
  if (parallel) log.printf("  parallelizing this calculation");
  domain.resize(2); domain[0]=d[0]; domain[1]=d[1]; nderivatives=nd;
}

void ActionWithValue::noAnalyticalDerivatives(){
  numericalDerivatives=true;
  warning("Numerical derivatives will be used as analytical derivatives are not available");
}

void ActionWithValue::gatherAllValues(){
  for(unsigned i=0;i<values.size();++i){
     // PARALLEL gather values[i]->value
     setValue( i, values[i]->value, 1.0 );
  }
}

void ActionWithValue::addValue( const std::string& name, const bool& ignorePeriod, const bool& hasDerivatives ){
  std::string thename; thename=getLabel() + "." + name;
  for(unsigned i=0;i<values.size();++i){
     if( values[i]->myname==thename ) assert(false);
  }
  
  unsigned nder;
  if( hasDerivatives ){ nder=nderivatives; } else{ nder=0; }
  if( ignorePeriod ){
      std::vector<double> fdomain(2,0.0);
      values.push_back(new Value(*this, thename, nder, fdomain) );
  } else {
      values.push_back(new Value(*this, thename, nder, domain) );
  } 
}

void ActionWithValue::clearInputForces(){
  for(unsigned i=0;i<values.size();i++){
     values[i]->hasForce=false; values[i]->inputForce=0.0;
  }
}

void ActionWithValue::clearDerivatives(){
  for(unsigned i=0;i<values.size();i++) values[i]->clearDerivatives();
}

Value* ActionWithValue::getValuePointer( const unsigned& i ){ 
  assert( i<values.size() );
  return values[i]; 
}

Value* ActionWithValue::getValuePointer( const std::string& name ){
  std::string thename=getLabel() + "." +  name;
  for(unsigned i=0;i<values.size();++i){
     if( values[i]->myname==thename ) return values[i];
  }
  error("there is no component with label " + name );
  return values[0];
}
