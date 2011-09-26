#ifndef __PLUMED_Value_h
#define __PLUMED_Value_h

#include <vector>
#include <string>
#include <cassert>


namespace PLMD{

class ActionWithValue;
class ActionWithArguments;

/// Class containing a value which can be addressed by PLMD::ActionWithArguments.
/// It also contains the derivative of this value with respect to
/// an arbitrary number of parameters.
/// Typically, an object of type PLMD::ActionWithValue will contain one or
/// more objects of type PLUMD::Value, one per component.
class Value{
friend class ActionWithValue;
friend class ActionWithArguments;
private:
/// The action that calculates this particular value
  ActionWithValue& action;
/// The name of this component
  std::string myname;
/// The positionVector
  double value;
/// Is the vector space we are mapping into periodic
  enum {unset,periodic,notperiodic} periodicity;
// The minimum and maximum for the periodicity in our new vector space
  double min,max;
/// The force on this quantity
  double inputForce;
/// The derivatives of the action
  std::vector<double> derivatives;
/// Is there a chain rule force on this quantity
  bool hasForce;
/// Do we have derivatives for this quantity
  bool deriv;
/// Clear all the derivatives
  void clearDerivatives();
/// Set the value of the quantity (this is only ever used ActionWithArguments::numericalDerivatives)
  void set( const double& f );
/// Get the value of a particular numbered derivative
  double getDerivative( const unsigned& j ) const ; 
public:
  Value( ActionWithValue&action, const std::string& name, const unsigned& nd, const std::vector<double>& domain );
/// Get the number of derivatives
  unsigned getNumberOfDerivatives() const ;
/// Calculate a difference 
  double difference( double c) const;
  double difference( double, double ) const;
/// Get the action that corresponds to this action
  ActionWithValue& getAction() const ;
/// Get the forces
  bool getForces( std::vector<double>& forces ) const ;
};

inline
unsigned Value::getNumberOfDerivatives() const { return derivatives.size(); }

inline
double Value::getDerivative( const unsigned& j ) const {
  assert( j<derivatives.size() );
  return derivatives[j]; 
}

inline
void Value::clearDerivatives(){
  derivatives.assign(derivatives.size(),0.0);
}

inline
double Value::difference(double c) const { 
   return difference( value,c ); 
}

inline
void Value::set(const double& f ){
   value=f;
} 

inline
bool Value::getForces( std::vector<double>& forces ) const {
   if( !hasForce ) return false;
   assert( derivatives.size()==forces.size() );
   for(unsigned i=0;i<derivatives.size();++i){ forces[i]=inputForce*derivatives[i]; }
   return true;
}

inline
ActionWithValue& Value::getAction() const {
  return action;
}

}
#endif

