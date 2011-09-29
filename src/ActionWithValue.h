#ifndef __PLUMED_ActionWithValue_h
#define __PLUMED_ActionWithValue_h

#include "Action.h"
#include "Value.h"
#include <vector>
#include <cassert>

namespace PLMD{

class ActionWithArguments;

/// Used to store the values calculated by a given action.
/// This object contains an array of PLMD::Value, one for each component.
/// It also stores all the derivatives of these values wrt the parameters
/// Parameters are other values (from other Action s) or atomic positions.
class ActionWithValue : public Action {
friend class ActionWithArguments;
private:
/// The number of derivatives the values in this container should have
  unsigned nderivatives;
/// The periodicities of the values in this particular container
  std::vector<double> domain;
/// Are we using numerical derivatives
  bool numericalDerivatives;
/// The array of values
  std::vector<Value*> values;
/// Return a pointer to one of the values in this action
/// These two routines are ONLY ever used by ActionWithArguments
  Value* getValuePointer( const unsigned& i );
  Value* getValuePointer( const std::string& name );
protected:
/// If, for whatever reasons, you do not have analytical derivatives implemented in your function this routine 
/// will ensure that numerical derivatives will be calculated at all times
  void noAnalyticalDerivatives();
/// Read all keywords relevant to ActionWithValue
  void readActionWithValue( const unsigned& nd, const std::vector<double>& domain );
/// Used to make sure that numerical derivatives does not numerically differentiate quantities that are not differentiable
  bool isMeaningfulToDifferentiate( const unsigned& n ) const ;
/// Add a value to the action
  void addValue( const std::string& name, const bool& ignorePeriod, const bool& hasDerivatives );
/// Accumulate the derivatives
  void addDerivative( const unsigned& n, const unsigned& nd, const double& val );
/// Get the forces acting on a particular value ( force*derivatives )
  bool getForces( const unsigned& n, std::vector<double>& forces ) const ;
/// Set the value (also applies chain rule by multiplying all derivatives by df)
  void setValue( const unsigned& n, const double& f, const double& df );
/// Set the named value (also applies chain rule by multiplying all derivatives by df)
  void setValue( const std::string& name, const double& f, const double& df );
/// Return the index for the value with name valname
  unsigned getValueNumberForLabel( const std::string& valname );
/// Get the number of Values
  unsigned getNumberOfValues() const;
/// Get the value of a particular numbered value
  double getValue( const unsigned i ) const;
public:
  ActionWithValue(const ActionOptions&ao);
  ~ActionWithValue();
/// Clear the forces on the values
  void clearInputForces();
/// Clear the derivatives of values wrt parameters
  void clearDerivatives();
/// Check if numerical derivatives should be used
  bool checkNumericalDerivatives() const;
};

inline
unsigned ActionWithValue::getValueNumberForLabel( const std::string& valname ){
  std::string thename = getLabel() + "." + valname;
  for(unsigned i=0;i<values.size();++i){
     if( values[i]->myname==thename ) return i;
  }
  assert(false);
  return 0;
}

inline
void ActionWithValue::setValue( const unsigned& n, const double& f, const double& df ){
  assert( n<values.size() );
  values[n]->value=f;
  for(unsigned i=0;i<values[n]->derivatives.size();++i){ values[n]->derivatives[i]*=df; } 
}

inline
void ActionWithValue::setValue( const std::string& name, const double& f, const double& df ){
  unsigned n=getValueNumberForLabel( name );
  values[n]->value=f;
  for(unsigned i=0;i<values[n]->derivatives.size();++i){ values[n]->derivatives[i]*=df; }
}

inline
bool ActionWithValue::isMeaningfulToDifferentiate( const unsigned& n ) const {
   return values[n]->deriv;
}

inline
bool ActionWithValue::checkNumericalDerivatives()const{
  return numericalDerivatives;
}

inline
unsigned ActionWithValue::getNumberOfValues() const {
  return values.size();
}

inline
double ActionWithValue::getValue(const unsigned i) const {
  return values[i]->value;
}

inline
void ActionWithValue::addDerivative( const unsigned& n, const unsigned& nd, const double& val ){
  assert( n<values.size() ); assert( nd<values[n]->derivatives.size() );
  values[n]->derivatives[nd]+=val; 
}

inline
bool ActionWithValue::getForces( const unsigned& n, std::vector<double>& forces ) const {
  assert( n<values.size() );
  return values[n]->getForces( forces );  
}

}

#endif
