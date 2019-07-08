/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2019 The plumed team
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

class OFile;
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
  friend class PlumedMain;
  friend class ActionWithValue;
  friend class ActionWithArguments;
private:
/// Is this value created by plumedmain
  bool created_in_plumedmain;
/// The action in which this quantity is calculated
  ActionWithValue* action;
/// Had the value been set
  bool value_set, reset;
/// The value of the normalization
  double norm;
/// The data for this value
  std::vector<double> data;
/// The value of the quantity
//  double value;
/// The force acting on this quantity
  std::vector<double> inputForces;
/// A flag telling us we have a force acting on this quantity
  bool hasForce;
/// The derivatives of the quantity stored in value
//  std::vector<double> derivatives;
  std::map<AtomNumber,Vector> gradients;
/// The name of this quantiy
  std::string name;
/// Does this quanity have derivatives
  bool hasDeriv;
/// Is this value a time series
  bool istimeseries;
/// What is the shape of the value (0 dimensional=scalar, 1 dimensional=vector, 2 dimensional=matrix)
  std::vector<unsigned> shape;
/// This is used by actions that always store data.  They cannot operate without storing all values
  bool alwaysstore;
/// Are we storing the data
  bool storedata, neverstore;
  std::vector<std::pair<std::string,int> > store_data_for;
/// Are we taking column sums
  bool columnsums;
/// Some variables for dealing with matrices
  bool symmetric;
/// Variables for storing data
  unsigned bufstart, streampos, matpos;
/// Store information on who is using information contained in this value
  std::map<std::string,std::vector<std::pair<int,int> > > userdata;
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
/// A constructor that is used throughout the code to setup the value poiters
  Value(ActionWithValue* av, const std::string& name, const bool withderiv,const std::vector<unsigned>&ss=std::vector<unsigned>());
/// Add information on who is using this action
  void interpretDataRequest( const std::string& uselab, unsigned& nargs, std::vector<Value*>& args, const std::string& values );
/// Set the shape of the Value
  void setShape( const std::vector<unsigned>&ss );
/// Set the value of the function
  void set(double);
/// Set the value of the stored data
  void set(const unsigned& n, const double& v );
/// Add something to the value of the function
  void add(double);
/// Add something to the stored data
  void add(const unsigned& n, const double& v );
/// Get the value of the function
  double get() const;
/// Get the value of a particular function
  double get( const unsigned& ival ) const ;
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
/// Get the size of this value
  unsigned getNumberOfValues( const std::string& actlab ) const ;
/// Set the number of derivatives
  void resizeDerivatives(int n);
/// Set all the derivatives to zero
  void clearDerivatives();
/// Add some derivative to the ith component of the derivatives array
  void addDerivative(unsigned i,double d);
/// Set the value of the ith component of the derivatives array
  void setDerivative(unsigned i, double d);
/// Get the derivative with respect to component n
  double getDerivative(const unsigned n) const;
/// get the derivative of a grid at a point n with resepct to argument j
  double getGridDerivative(const unsigned& n, const unsigned& j ) const ;
/// Clear the input force on the variable
  void clearInputForce();
/// Add some force on this value
  void  addForce(const unsigned& iforce, double f);
/// Get the value of the force on this colvar
  double getForce( const unsigned& iforce ) const ;
/// Check if forces have been added on this value
  bool forcesWereAdded() const ;
/// Apply the forces to the derivatives using the chain rule (if there are no forces this routine returns false)
  bool applyForce( std::vector<double>& forces ) const ;
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
/// Build the store of data
  void buildDataStore( const std::string& actlabel );
  void makeTimeSeries();
  bool isTimeSeries() const ;
  void alwaysStoreValues();
  void neverStoreValues();
  void buildColumnSums();
///
  unsigned getRank() const ;
///
  double getNorm() const ;
///
  void setNorm( const double& nn );
///
  const std::vector<unsigned>& getShape() const ;
///
  unsigned getSize() const ;
///
  unsigned getPositionInStream() const ;
///
  unsigned getPositionInMatrixStash() const ;
///
  std::size_t getIndex(const std::vector<unsigned> & indices) const ;
///
  void convertIndexToindices(const std::size_t& index, std::vector<unsigned>& indices ) const ;
///
  void print( const std::string& alabel, OFile& ofile ) const ;
///
  bool usingAllVals( const std::string& alabel ) const ;
///
  double getRequiredValue(  const std::string& alabel, const unsigned& num  ) const ;
///
  void getRequiredValue( const std::string& alabel, const unsigned& num, std::vector<double>& args ) const ;
///
  std::string getOutputDescription( const std::string& alabel ) const ;
///
  std::string getOutputDescription( const std::string& alabel, const unsigned& i ) const ;
///
  void setSymmetric( const bool& sym );
///
  bool isSymmetric() const ;
};

inline
void Value::applyPeriodicity( const unsigned& ival ) {
  if(periodicity==periodic) {
    data[ival]=min+difference(min,data[ival]);
    if(data[ival]<min)data[ival]+=max_minus_min;
  }
}

inline
void Value::set(double v) {
  plumed_dbg_assert( shape.size()==0 );
  value_set=true;
  data[0]=v;
  applyPeriodicity(0);
}

inline
void Value::add(const unsigned& n, const double& v ) {
  value_set=true; data[n]+=v; applyPeriodicity(n);
}

inline
void Value::add(double v) {
  plumed_dbg_assert( shape.size()==0 );
  value_set=true;
  data[0]+=v;
  applyPeriodicity(0);
}

inline
double Value::get()const {
  plumed_dbg_assert( shape.size()==0 );
  if( norm>epsilon ) return data[0] / norm;
  return 0.0;
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
  plumed_massert(hasDeriv,"the derivatives array for value " + name + " has zero size" );
  if( shape.size()>0 ) return shape.size();
  return data.size()-1;
}

inline
double Value::getDerivative(const unsigned n) const {
  plumed_dbg_massert(shape.size()==0 && n<data.size()-1,"you are asking for a derivative that is out of bounds");
  return data[1+n] / norm;
}

inline
void Value::resizeDerivatives(int n) {
  if( shape.size()>0 ) return;
  if(hasDeriv) data.resize(1+n);
}

inline
void Value::addDerivative(unsigned i,double d) {
  plumed_dbg_massert( shape.size()==0 && i<data.size()-1,"derivative is out of bounds");
  data[1+i]+=d;
}

inline
void Value::setDerivative(unsigned i, double d) {
  plumed_dbg_massert( shape.size()==0 && i<data.size()-1,"derivative is out of bounds");
  data[1+i]=d;
}

inline
void Value::clearInputForce() {
  hasForce=false;
  std::fill(inputForces.begin(),inputForces.end(),0);
}

inline
void Value::clearDerivatives() {
  value_set=false; if( data.size()>1 ) std::fill(data.begin()+1, data.end(), 0);
}

inline
void Value::addForce(const unsigned& iforce, double f) {
  plumed_dbg_assert( iforce<inputForces.size() );
  hasForce=true;
  inputForces[iforce]+=f;
}

inline
double Value::getForce( const unsigned& iforce ) const {
  plumed_dbg_assert( iforce<inputForces.size() );
  return inputForces[iforce];
}

inline
bool Value::forcesWereAdded() const {
  return hasForce;
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
  } else plumed_merror("periodicity in " + name + " should be set to compute differences");
}

inline
double Value::bringBackInPbc(double d1)const {
  return min+max_minus_min/2.+difference(min+max_minus_min/2., d1);
}

inline
unsigned Value::getRank() const {
  return shape.size();
}

inline
double Value::getNorm() const {
  return norm;
}

inline
void Value::setNorm( const double& nn ) {
  norm = nn;
}


inline
std::size_t Value::getIndex(const std::vector<unsigned> & indices) const {
  plumed_dbg_assert( indices.size()==shape.size() );
#ifndef DNDEBUG
  for(unsigned i=0; i<indices.size(); ++i) plumed_dbg_assert( indices[i]<shape[i] );
#endif
  std::size_t index=indices[shape.size()-1];
  for(unsigned i=shape.size()-1; i>0 ; --i ) index=index*shape[i-1]+indices[i-1];
  return index;
}

inline
void Value::convertIndexToindices(const std::size_t& index, std::vector<unsigned>& indices ) const {
  std::size_t kk=index; indices[0]=index%shape[0];
  for(unsigned i=1; i<shape.size()-1; ++i) {
    kk=(kk-indices[i-1])/shape[i-1];
    indices[i]=kk%shape[i];
  }
  if(shape.size()>=2) indices[shape.size()-1]=(kk-indices[shape.size()-2])/shape[shape.size()-2];
}

inline
double Value::getMaxMinusMin()const {
  plumed_dbg_assert( periodicity==periodic );
  return max_minus_min;
}

inline
bool Value::isTimeSeries() const {
  return istimeseries;
}

}

#endif

