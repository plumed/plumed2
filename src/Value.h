#ifndef __PLUMED_Value_h
#define __PLUMED_Value_h

#include <vector>
#include <string>
#include <map>
#include "PlumedException.h"
#include "Tools.h"
#include "AtomNumber.h"
#include "Vector.h"

namespace PLMD{

class ActionWithValue;

//+DEVELDOC TOOLBOX Value
/**
A class for holding the value of a function together with its derivatives.
*/
//+ENDDEVELDOC

/// Typically, an  object of type PLMD::ActionWithValue will contain one 
/// object of type PLUMD::Value that will be named after the label.  If the 
/// PLMD::ActionWithValue is part of a class that calculates multiple components 
/// then the class will contain multiple that will be called label.component-name
/// This class is used to pass information between different PLMD::Action 
/// objects.  However, if you find a use for a tempory PLMD::Value in some method
/// you are implementing please feel free to use it.
class Value{
friend class ActionWithValue;
private:
/// The action in which this quantity is calculated
  ActionWithValue* action;
/// Had the value been set
  bool value_set;
/// The value of the quantity
  double value;
/// The force acting on this quantity
  double inputForce;
/// A flag telling us we have a force acting on this quantity
  bool hasForce;
/// The derivatives of the quantity stored in value
  std::vector<double> derivatives;
  std::map<AtomNumber,Vector> gradients;
/// The name of this quantiy
  std::string name;
/// Does this quanity have derivatives
  bool hasDeriv;
/// Is this quantity periodic
  enum {unset,periodic,notperiodic} periodicity;
/// Various quantities that describe the domain of this value
  double min,max;
  double max_minus_min;
  double inv_max_minus_min;
/// Complete the setup of the periodicity
  void setupPeriodicity();
public:
/// A constructor that can be used to make Vectors of values
  Value() : action(NULL), hasDeriv(true) {} ;
/// A constructor that is used throughout the code to setup the value poiters
  Value(ActionWithValue* av, const std::string& name, const bool withderiv);
/// Set the value of the function 
  void set(double);
/// Get the value of the function
  double get() const;
  void setPeriodicity(bool);
  void setDomain(double,double);
/// Check if the value is periodic
  bool isPeriodic() const;
/// Get the domain of the quantity
  void getDomain(double&,double&) const;
/// Get the name of the quantity
  const std::string& getName() const;
/// Check whether or not this particular quantity has derivatives
  bool hasDerivatives()const;
/// Get the number of derivatives that this particular value has
  unsigned getNumberOfDerivatives() const; 
/// Set the number of derivatives
  void resizeDerivatives(int n);
/// Set all the derivatives to zero
  void clearDerivatives();
/// Add some derivative to the ith component of the derivatives array
  void addDerivative(int i,double d);
/// Get the derivative with respect to component n
  double getDerivative(const unsigned n) const;
/// Clear the input force on the variable
  void clearInputForce();
/// Add some force on this value
  void  addForce(double f);
/// Get the value of the force on this colvar
  double getForce() const ;
/// Apply the forces to the derivatives using the chain rule (if there are no forces this routine returns false
  bool applyForce( std::vector<double>& forces ) const ;

/// Calculate the difference between the instantaneous value of the function and some other point
  double difference(double)const;
/// Calculate the difference between two values of this function
  double difference(double,double)const;
/// This sets up the gradients
  void setGradients();
  static double projection(const Value&,const Value&);
};

inline
void Value::set(double v){
  value=v;
}

inline
double Value::get()const{
  return value;
}

inline
const std::string& Value::getName()const{
  return name;
}

inline
unsigned Value::getNumberOfDerivatives() const {
  plumed_massert(derivatives.size()!=0,"the derivatives array for this value has zero size");
  return derivatives.size();
}

inline
double Value::getDerivative(const unsigned n) const {
  plumed_massert(n<derivatives.size(),"you are asking for a derivative that is out of bounds");
  return derivatives[n];
}

inline
bool Value::hasDerivatives() const {
  return (derivatives.size()!=0);
}

inline
void Value::resizeDerivatives(int n){
  plumed_massert(hasDeriv,"cannot resize derivatives in values that have not got derivatives"); 
  derivatives.resize(n);
}

inline
void Value::addDerivative(int i,double d){
  plumed_massert(i<derivatives.size(),"derivative is out of bounds");
  derivatives[i]=d;
}

inline
void Value::clearInputForce(){
  hasForce=false;
  inputForce=0.0;
}

inline
void Value::clearDerivatives(){
  derivatives.assign(derivatives.size(),0.0);
}

inline
void Value::addForce(double f){
  plumed_massert(hasDerivatives(),"forces can only be added to values with derivatives");
  hasForce=true;
  inputForce+=f;
}

inline
double Value::getForce() const {
  return inputForce;
}

inline
double Value::difference(double d1,double d2)const{
  if(periodicity==notperiodic){
    return d2-d1;
  }else if(periodicity==periodic){
    double s=(d2-d1)*inv_max_minus_min;
    s=Tools::pbc(s);
    return s*max_minus_min;
  } else plumed_merror("periodicity should be set to compute differences");
}

inline
double Value::difference(double d)const{
  return difference(get(),d);
}

}

#endif

