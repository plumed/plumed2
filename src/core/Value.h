/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
  friend class DomainDecomposition;
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
  std::vector<double> data;
/// The force acting on this quantity
  std::vector<double> inputForce;
/// A flag telling us we have a force acting on this quantity
  bool hasForce;
/// This flag is used if the value is a constant
  bool constant;
/// The derivatives of the quantity stored in value
  std::map<AtomNumber,Vector> gradients;
/// The name of this quantiy
  std::string name;
/// Are we storing the data for this value if it is vector or matrix
  bool storedata;
/// What is the shape of the value (0 dimensional=scalar, n dimensional with derivatives=grid, 1 dimensional no derivatives=vector, 2 dimensional no derivatives=matrix)
  std::vector<unsigned> shape;
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
  void applyPeriodicity( const unsigned& ival );
public:
/// A constructor that can be used to make Vectors of values
  Value();
/// A constructor that can be used to make Vectors of named values
  explicit Value(const std::string& name);
/// A constructor that is used throughout the code to setup the value poiters
  Value(ActionWithValue* av, const std::string& name, const bool withderiv,const std::vector<unsigned>&ss=std::vector<unsigned>());
/// Set the shape of the Value
  void setShape( const std::vector<unsigned>&ss );
/// Set the value of the function
  void set(double);
/// Set the value of the stored data
  void set(const unsigned& n, const double& v );
/// Add something to the value of the function
  void add(double);
/// Get the value of the function
  double get( const unsigned& ival=0 ) const;
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
  void addForce(double f);
/// Add some force on the ival th component of this value
  void addForce( const unsigned& ival, double f );
/// Get the value of the force on this colvar
  double getForce( const unsigned& ival=0 ) const ;
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
/// Get the rank of the object that is contained in this value
  unsigned getRank() const ;
/// Get the shape of the object that is contained in this value
  const std::vector<unsigned>& getShape() const ;
/// This turns on storing of vectors/matrices
  void buildDataStore();
/// Get the total number of scalars that are stored here
  unsigned getNumberOfValues() const ;
/// Get the number of threads to use when assigning this value
  unsigned getGoodNumThreads( const unsigned& j, const unsigned& k ) const ;
/// These are used for passing around the data in this value when we are doing replica exchange
  void writeBinary(std::ostream&o) const ;
  void readBinary(std::istream&i);
/// These are used for making constant values
  bool isConstant() const ;
  void setConstant();
/// Check if forces have been added on this value
  bool forcesWereAdded() const ;
};

void copy( const Value& val1, Value& val2 );
void copy( const Value& val1, Value* val2 );
void add( const Value& val1, Value* valout );

inline
void Value::applyPeriodicity(const unsigned& ival) {
  if(periodicity==periodic) {
    data[ival]=min+difference(min,data[ival]);
    if(data[ival]<min)data[ival]+=max_minus_min;
  }
}

inline
void product( const Value& val1, const Value& val2, Value& valout ) {
  plumed_assert( val1.getNumberOfDerivatives()==val2.getNumberOfDerivatives() );
  if( valout.getNumberOfDerivatives()!=val1.getNumberOfDerivatives() ) valout.resizeDerivatives( val1.getNumberOfDerivatives() );
  valout.value_set=false;
  valout.clearDerivatives();
  double u=val1.get();
  double v=val2.get();
  for(unsigned i=0; i<val1.getNumberOfDerivatives(); ++i) {
    valout.addDerivative(i, u*val2.getDerivative(i) + v*val1.getDerivative(i) );
  }
  valout.set( u*v );
}

inline
void quotient( const Value& val1, const Value& val2, Value* valout ) {
  plumed_assert( val1.getNumberOfDerivatives()==val2.getNumberOfDerivatives());
  if( valout->getNumberOfDerivatives()!=val1.getNumberOfDerivatives() ) valout->resizeDerivatives( val1.getNumberOfDerivatives() );
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
  data[0]=v;
  applyPeriodicity(0);
}

inline
void Value::add(double v) {
  value_set=true;
  data[0]+=v;
  applyPeriodicity(0);
}

inline
double Value::get( const unsigned& ival )const {
  return data[ival];
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
  if( shape.size()>0 ) return shape.size();
  return data.size() - 1;
}

inline
double Value::getDerivative(const unsigned n) const {
  plumed_dbg_massert(n<getNumberOfDerivatives(),"you are asking for a derivative that is out of bounds");
  return data[1+n];
}

inline
bool Value::hasDerivatives() const {
  return hasDeriv;
}

inline
void Value::resizeDerivatives(int n) {
  if( shape.size()>0 ) return;
  if(hasDeriv) data.resize(1+n);
}

inline
void Value::addDerivative(unsigned i,double d) {
  plumed_dbg_massert(i<getNumberOfDerivatives(),"derivative is out of bounds");
  data[1+i]+=d;
}

inline
void Value::setDerivative(unsigned i, double d) {
  plumed_dbg_massert(i<getNumberOfDerivatives(),"derivative is out of bounds");
  data[1+i]=d;
}

inline
void Value::chainRule(double df) {
  for(unsigned i=0; i<getNumberOfDerivatives(); ++i) data[1+i]*=df;
}

inline
void Value::clearInputForce() {
  hasForce=false;
  std::fill(inputForce.begin(),inputForce.end(),0);
}

inline
void Value::clearDerivatives() {
  if( constant ) return;
  value_set=false;
  if( data.size()>1 ) std::fill(data.begin()+1, data.end(), 0);
}

inline
void Value::addForce(double f) {
  hasForce=true;
  inputForce[0]+=f;
}

inline
void Value::addForce(const unsigned& ival, double f) {
  plumed_dbg_massert(ival<inputForce.size(),"too few components in value to add force");
  hasForce=true;
  inputForce[ival]+=f;
}

inline
bool Value::forcesWereAdded() const {
  return hasForce;
}

inline
double Value::getForce( const unsigned& ival ) const {
  plumed_dbg_assert( ival<inputForce.size() );
  return inputForce[ival];
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

inline
unsigned Value::getRank() const {
  return shape.size();
}

inline
const std::vector<unsigned>& Value::getShape() const {
  return shape;
}

inline
unsigned Value::getNumberOfValues() const {
  unsigned size=1; for(unsigned i=0; i<shape.size(); ++i) size *= shape[i];
  return size;
}

inline
bool Value::isConstant() const {
  return constant;
}

}

#endif

