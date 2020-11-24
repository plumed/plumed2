/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2020 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#ifndef __PLUMED_core_Value_h
#define __PLUMED_core_Value_h

#include <vector>
#include <string>
#include <map>
#include "tools/Exception.h"
#include "tools/Tools.h"
#include "tools/AtomNumber.h"
#include "tools/Vector.h"

namespace PLMD {

class ActionWithValue;

/// \ingroup TOOLBOX
/// A class for holding the value of a function together with its derivatives.
/// Typically, an  object of type PLMD::ActionWithValue will contain one
/// object of type PLUMD::Value that will be named after the label.  If the
/// PLMD::ActionWithValue is part of a class that calculates multiple components
/// then the class will contain multiple that will be called label.component-name
/// This class is used to pass information between different PLMD::Action
/// objects.  However, if you find a use for a tempory PLMD::Value in some method
/// you are implementing please feel free to use it.
class Value {
  friend class ActionWithValue;
/// This copies the contents of a value into a second value (just the derivatives and value)
  friend void copy( const Value& val1, Value& val2 );
/// This copies the contents of a value into a second value (but second value is a pointer)
  friend void copy( const Value& val, Value* val2 );
/// This adds some derivatives onto the value
  friend void add( const Value& val1, Value* valout );
/// This calculates val1*val2 and sorts out the derivatives
  friend void product( const Value& val1, const Value& val2, Value& valout );
/// This calculates va1/val2 and sorts out the derivatives
  friend void quotient( const Value& val1, const Value& val2, Value* valout );
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
  std::string str_min, str_max;
  double min,max;
  double max_minus_min;
  double inv_max_minus_min;
/// Complete the setup of the periodicity
  void setupPeriodicity();
// bring value within PBCs
  void applyPeriodicity();
public:
/// A constructor that can be used to make Vectors of values
  Value();
/// A constructor that can be used to make Vectors of named values
  explicit Value(const std::string& name);
/// A constructor that is used throughout the code to setup the value poiters
  Value(ActionWithValue* av, const std::string& name, const bool withderiv);
/// Set the value of the function
  void set(double);
/// Add something to the value of the function
  void add(double);
/// Get the value of the function
  double get() const;
/// Find out if the value has been set
  bool valueHasBeenSet() const;
/// Check if the value is periodic
  bool isPeriodic() const;
/// Set the function not periodic
  void setNotPeriodic();
/// Set the domain of the function
  void setDomain(const std::string&, const std::string&);
/// Get the domain of the quantity
  void getDomain(std::string&,std::string&) const;
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
  void addDerivative(unsigned i,double d);
/// Set the value of the ith component of the derivatives array
  void setDerivative(unsigned i, double d);
/// Apply the chain rule to the derivatives
  void chainRule(double df);
/// Get the derivative with respect to component n
  double getDerivative(const unsigned n) const;
/// Clear the input force on the variable
  void clearInputForce();
/// Add some force on this value
  void  addForce(double f);
/// Get the value of the force on this colvar
  double getForce() const ;
/// Apply the forces to the derivatives using the chain rule (if there are no forces this routine returns false)
  bool applyForce( std::vector<double>& forces ) const ;
/// Calculate the difference between the instantaneous value of the function and some other point: other_point-inst_val
  double difference(double)const;
/// Calculate the difference between two values of this function: d2 -d1
  double difference(double d1,double d2)const;
/// This returns the pointer to the action where this value is calculated
  ActionWithValue* getPntrToAction();
/// Bring back one value into the correct pbc if needed, else give back the value
  double bringBackInPbc(double d1)const;
/// Get the difference between max and minimum of domain
  double getMaxMinusMin()const;
/// This sets up the gradients
  void setGradients();
  static double projection(const Value&,const Value&);
};

void copy( const Value& val1, Value& val2 );
void copy( const Value& val1, Value* val2 );
void add( const Value& val1, Value* valout );

inline
void Value::applyPeriodicity() {
  if(periodicity==periodic) {
    value=min+difference(min,value);
    if(value<min)value+=max_minus_min;
  }
}

inline
void product( const Value& val1, const Value& val2, Value& valout ) {
  plumed_assert( val1.derivatives.size()==val2.derivatives.size() );
  if( valout.derivatives.size()!=val1.derivatives.size() ) valout.resizeDerivatives( val1.derivatives.size() );
  valout.value_set=false;
  valout.clearDerivatives();
  double u=val1.value;
  double v=val2.value;
  for(unsigned i=0; i<val1.derivatives.size(); ++i) {
    valout.addDerivative(i, u*val2.derivatives[i] + v*val1.derivatives[i] );
  }
  valout.set( u*v );
}

inline
void quotient( const Value& val1, const Value& val2, Value* valout ) {
  plumed_assert( val1.derivatives.size()==val2.derivatives.size() );
  if( valout->derivatives.size()!=val1.derivatives.size() ) valout->resizeDerivatives( val1.derivatives.size() );
  valout->value_set=false;
  valout->clearDerivatives();
  double u=val1.get();
  double v=val2.get();
  for(unsigned i=0; i<val1.getNumberOfDerivatives(); ++i) {
    valout->addDerivative(i, v*val1.getDerivative(i) - u*val2.getDerivative(i) );
  }
  valout->chainRule( 1/(v*v) ); valout->set( u / v );
}

inline
void Value::set(double v) {
  value_set=true;
  value=v;
  applyPeriodicity();
}

inline
void Value::add(double v) {
  value_set=true;
  value+=v;
  applyPeriodicity();
}

inline
double Value::get()const {
  return value;
}

inline
bool Value::valueHasBeenSet() const {
  return value_set;
}

inline
const std::string& Value::getName()const {
  return name;
}

inline
unsigned Value::getNumberOfDerivatives() const {
  plumed_massert(hasDeriv,"the derivatives array for this value has zero size");
  return derivatives.size();
}

inline
double Value::getDerivative(const unsigned n) const {
  plumed_dbg_massert(n<derivatives.size(),"you are asking for a derivative that is out of bounds");
  return derivatives[n];
}

inline
bool Value::hasDerivatives() const {
  return hasDeriv;
}

inline
void Value::resizeDerivatives(int n) {
  if(hasDeriv) derivatives.resize(n);
}

inline
void Value::addDerivative(unsigned i,double d) {
  plumed_dbg_massert(i<derivatives.size(),"derivative is out of bounds");
  derivatives[i]+=d;
}

inline
void Value::setDerivative(unsigned i, double d) {
  plumed_dbg_massert(i<derivatives.size(),"derivative is out of bounds");
  derivatives[i]=d;
}

inline
void Value::chainRule(double df) {
  for(unsigned i=0; i<derivatives.size(); ++i) derivatives[i]*=df;
}

inline
void Value::clearInputForce() {
  hasForce=false;
  inputForce=0.0;
}

inline
void Value::clearDerivatives() {
  value_set=false;
  std::fill(derivatives.begin(), derivatives.end(), 0);
}

inline
void Value::addForce(double f) {
  plumed_dbg_massert(hasDerivatives(),"forces can only be added to values with derivatives");
  hasForce=true;
  inputForce+=f;
}

inline
double Value::getForce() const {
  return inputForce;
}
/// d2-d1
inline
double Value::difference(double d1,double d2)const {
  if(periodicity==notperiodic) {
    return d2-d1;
  } else if(periodicity==periodic) {
    double s=(d2-d1)*inv_max_minus_min;
    // remember: pbc brings the difference in a range of -0.5:0.5
    s=Tools::pbc(s);
    return s*max_minus_min;
  } else plumed_merror("periodicity should be set to compute differences");
}

inline
double Value::bringBackInPbc(double d1)const {
  return min+max_minus_min/2.+difference(min+max_minus_min/2., d1);
}

inline
double Value::difference(double d)const {
  return difference(get(),d);
}

inline
double Value::getMaxMinusMin()const {
  plumed_dbg_assert( periodicity==periodic );
  return max_minus_min;
}

}

#endif

