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
#include "tools/View.h"

namespace PLMD {

class OFile;
class Communicator;
class ActionWithValue;
class ActionAtomistic;

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
  friend struct ArgumentsBookkeeping;
  friend class ActionWithValue;
  friend class ActionWithVector;
  friend class ActionWithMatrix;
  friend class ActionAtomistic;
  friend class ActionWithArguments;
  friend class ActionWithVirtualAtom;
  friend class DomainDecomposition;
  template<typename T>
  friend class DataPassingObjectTyped;
private:
/// The action in which this quantity is calculated
  ActionWithValue* action=nullptr;
/// Had the value been set
  bool value_set=false;
/// The value of the quantity
  std::vector<double> data;
/// The force acting on this quantity
  std::vector<double> inputForce;
/// A flag telling us we have a force acting on this quantity
  bool hasForce=false;
/// The way this value is used in the code
/// normal = regular value that is determined during calculate
/// constant = constnt value that is determined during startup and that doesn't change during simulation
/// average = value that is averaged/collected over multiple steps of trajectory
/// calcFromAverage = value that is calculated from an average value
  enum {normal,constant,average,calcFromAverage} valtype=normal;
/// This is used by ActionWithValue to set the valtype
  void setValType( const std::string& vtype );
/// This is used by ActionWithValue to determine if we need to calculate on update
  bool calculateOnUpdate() const ;
/// The derivatives of the quantity stored in value
  std::map<AtomNumber,Vector> gradients;
/// The name of this quantiy
  std::string name;
/// What is the shape of the value (0 dimensional=scalar, n dimensional with derivatives=grid, 1 dimensional no derivatives=vector, 2 dimensional no derivatives=matrix)
  std::vector<std::size_t> shape;
/// Does this quanity have derivatives
  bool hasDeriv=true;
/// Variables for storing data
  unsigned bufstart=0;
  unsigned ngrid_der=0;
  std::size_t ncols=0;
/// If we are storing a matrix is it symmetric?
  bool symmetric=false;
/// This is a bookeeping array that holds the non-zero elements of the "sparse" matrix
  std::vector<unsigned> matrix_bookeeping;
/// Is this quantity periodic
  enum {unset,periodic,notperiodic} periodicity=unset;
/// Various quantities that describe the domain of this value
  std::string str_min;
  std::string str_max;
  double min=0.0;
  double max=0.0;
  double max_minus_min=0.0;
  double inv_max_minus_min=0.0;
/// Is the derivative of this quantity zero when the value is zero
  bool derivativeIsZeroWhenValueIsZero=false;
/// Complete the setup of the periodicity
  void setupPeriodicity();
// bring value within PBCs
  void applyPeriodicity( const unsigned& ival );
public:
/// A constructor that can be used to make Vectors of values
  Value();
/// A constructor that can be used to make Vectors of named values
  explicit Value(const std::string& valname);
/// A constructor that is used throughout the code to setup the value poiters
  Value(ActionWithValue* av, const std::string& valname, const bool withderiv,const std::vector<std::size_t>&ss=std::vector<std::size_t>());
/// Set the shape of the Value
  void setShape( const std::vector<std::size_t>&ss );
/// Set the value of the function
  void set(double);
/// Set the value of the stored data
  void set(const std::size_t& n, const double& v );
/// Add something to the value of the function
  void add(double);
/// Add something to the ith element of the data array
  void add(const std::size_t& n, const double& v );
/// Get the location of this element of in the store
  std::size_t getIndexInStore( const std::size_t& ival ) const ;
/// Get the value of the function
  double get( const std::size_t ival=0, const bool trueind=true ) const;
/// A variant of get() for checking that at least one of the values on the row is !=0
  bool checkValueIsActiveForMMul(std::size_t task) const;
/// A variant of get() for assigning data to an external view (assuems trueind=false), returns the number of arguments assigned
  std::size_t assignValues(View<double> target);
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
  void clearDerivatives( const bool force=false );
/// Add some derivative to the ith component of the derivatives array
  void addDerivative(unsigned i,double d);
/// Set the value of the ith component of the derivatives array
  void setDerivative(unsigned i, double d);
/// Get the derivative with respect to component n
  double getDerivative(const unsigned n) const;
/// Clear the input force on the variable
  void clearInputForce();
/// Special method for clearing forces on variables used by DataPassingObject
  void clearInputForce( const std::vector<AtomNumber>& index );
/// Set hasForce equal to true
  void addForce();
/// Add some force on this value
  void addForce(double f);
/// Add some force on the ival th component of this value
  void addForce( const std::size_t& ival, double f, const bool trueind=true );
  ///Add forces from a vector, imples trueInd=false and retunrs the number of forces assigned
  std::size_t addForces(View<const double> f);
/// Get the value of the force on this colvar
  double getForce( const std::size_t& ival=0 ) const ;
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
  void setGradients( ActionAtomistic* aa, unsigned& start );
/// This passes gradients from one action to another
  void passGradients( const double& der, std::map<AtomNumber,Vector>& g ) const ;
  static double projection(const Value&,const Value&);
/// Get the rank of the object that is contained in this value
  unsigned getRank() const ;
/// Get the shape of the object that is contained in this value
  const std::vector<std::size_t>& getShape() const ;
/// Reshape the storage for sparse matrices
  void reshapeMatrixStore( const unsigned& n );
/// Copy the matrix bookeeping stuff
  void copyBookeepingArrayFromArgument( Value* myarg );
/// Set the symmetric flag equal true for this matrix
  void setSymmetric( const bool& sym );
/// Get the total number of scalars that are stored here
  std::size_t getNumberOfValues() const ;
/// Get the number of values that are actually stored here once sparse matrices are taken into account
  std::size_t getNumberOfStoredValues() const ;
/// Get the number of threads to use when assigning this value
  unsigned getGoodNumThreads( const unsigned& j, const unsigned& k ) const ;
/// These are used for passing around the data in this value when we are doing replica exchange
  void writeBinary(std::ostream&o) const ;
  void readBinary(std::istream&i);
/// These are used for making constant values
  bool isConstant() const ;
  void setConstant();
  void reshapeConstantValue( const std::vector<std::size_t>& sh );
/// Check if forces have been added on this value
  bool forcesWereAdded() const ;
/// Set a bool that tells us if the derivative is zero when the value is zero true
  void setDerivativeIsZeroWhenValueIsZero();
/// Return a bool that tells us if the derivative is zero when the value is zero
  bool isDerivativeZeroWhenValueIsZero() const ;
/// Convert the input index to its corresponding indices
  void convertIndexToindices(const std::size_t& index, std::vector<unsigned>& indices ) const ;
/// Print out all the values in this Value
  void print( OFile& ofile ) const ;
/// Print out all the forces in this Value
  void printForce( OFile& ofile ) const ;
/// Are we to ignore the stored value
  bool ignoreStoredValue(const std::string& n) const ;
/// Set a matrix element to be non zero
  void setMatrixBookeepingElement( const unsigned& i, const unsigned& n );
///
  unsigned getRowLength( const std::size_t& irow ) const ;
///
  unsigned getRowIndex( std::size_t irow, std::size_t jind ) const ;
///
  void setRowIndices( const std::size_t& irow, const std::vector<std::size_t>& ind );
///
  std::size_t getNumberOfColumns() const ;
///
  bool isSymmetric() const ;
/// Retrieve the non-zero edges in a matrix
  void retrieveEdgeList( unsigned& nedge, std::vector<std::pair<unsigned,unsigned> >& active, std::vector<double>& elems );
/// Get the number of derivatives that the grid has
  unsigned getNumberOfGridDerivatives() const ;
/// get the derivative of a grid at a point n with resepct to argument j
  double getGridDerivative(const unsigned& n, const unsigned& j ) const ;
/// Add the derivatives of the grid to the corner
  void addGridDerivatives( const unsigned& n, const unsigned& j, const double& val );
///
  void setGridDerivatives( const unsigned& n, const unsigned& j, const double& val );
/// Add another value to the end of the data vector held by this value.  This is used in COLLECT
  void push_back( const double& val );
/// Get the type of value that is stored here
  std::string getValueType() const ;
/// Check if all the elements in the value have the same value
  bool allElementsEqual() const ;
};

inline
void Value::applyPeriodicity(const unsigned& ival) {
  if(periodicity==periodic) {
    data[ival]=min+difference(min,data[ival]);
    if(data[ival]<min) {
      data[ival]+=max_minus_min;
    }
  }
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
void Value::add(const std::size_t& n, const double& v ) {
  value_set=true;
  data[n]+=v;
  applyPeriodicity(n);
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
  if( shape.size()>0 ) {
    return shape.size();
  }
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
  if( shape.size()>0 ) {
    return;
  }
  if(hasDeriv) {
    data.resize(1+n);
  }
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
void Value::clearInputForce() {
  if( !hasForce ) {
    return;
  }
  hasForce=false;
  std::fill(inputForce.begin(),inputForce.end(),0);
}

inline
void Value::clearInputForce( const std::vector<AtomNumber>& index ) {
  if( !hasForce ) {
    return;
  }
  hasForce=false;
  for(const auto & p : index) {
    inputForce[p.index()]=0;
  }
}

inline
void Value::clearDerivatives( const bool force ) {
  if( !force && (valtype==constant || valtype==average) ) {
    return;
  }

  value_set=false;
  if( shape.size()>0 ) {
    std::fill(data.begin(), data.end(), 0);
  } else if( data.size()>1 ) {
    std::fill(data.begin()+1, data.end(), 0);
  }
}

inline
void Value::addForce() {
  hasForce=true;
}

inline
void Value::addForce(double f) {
  hasForce=true;
  inputForce[0]+=f;
}

inline
bool Value::forcesWereAdded() const {
  return hasForce;
}

inline
double Value::getForce( const std::size_t& ival ) const {
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
  } else {
    plumed_merror("periodicity should be set to compute differences");
  }
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
  if( valtype==constant && shape.size()==1 && shape[0]==1 ) {
    return 0;
  }
  return shape.size();
}

inline
const std::vector<std::size_t>& Value::getShape() const {
  return shape;
}

inline
std::size_t Value::getNumberOfValues() const {
  std::size_t size=1;
  for(unsigned i=0; i<shape.size(); ++i) {
    size *= shape[i];
  }
  return size;
}

inline
std::size_t Value::getNumberOfStoredValues() const {
  if( getRank()==2 && !hasDeriv ) {
    return shape[0]*ncols;
  }
  return getNumberOfValues();
}

inline
bool Value::isConstant() const {
  return valtype==constant;
}

inline
void Value::setDerivativeIsZeroWhenValueIsZero() {
  derivativeIsZeroWhenValueIsZero=true;
}

inline
bool Value::isDerivativeZeroWhenValueIsZero() const {
  return derivativeIsZeroWhenValueIsZero;
}

inline
void Value::setMatrixBookeepingElement( const unsigned& i, const unsigned& n ) {
  plumed_dbg_assert( i<matrix_bookeeping.size() );
  matrix_bookeeping[i]=n;
}

inline
unsigned Value::getRowLength( const std::size_t& irow ) const {
  if( matrix_bookeeping.size()==0 ) {
    return 0;
  }
  plumed_dbg_assert( (1+ncols)*irow<matrix_bookeeping.size() );
  return matrix_bookeeping[(1+ncols)*irow];
}

inline
unsigned Value::getRowIndex(const std::size_t irow, const std::size_t jind ) const {
  plumed_dbg_massert( (1+ncols)*irow+1+jind<matrix_bookeeping.size() && jind<matrix_bookeeping[(1+ncols)*irow], "failing in value " + name );
  return matrix_bookeeping[(1+ncols)*irow+1+jind];
}

inline
void Value::setRowIndices( const std::size_t& irow, const std::vector<std::size_t>& ind ) {
  plumed_dbg_massert( (1+ncols)*irow+1+ind.size()<=matrix_bookeeping.size(), "problem in " + name );
  std::size_t istart = (1+ncols)*irow;
  matrix_bookeeping[istart] = ind.size();
  ++istart;
  for(unsigned i=0; i<ind.size(); ++i) {
    matrix_bookeeping[istart] = ind[i];
    ++istart;
  }
}

inline
std::size_t Value::getNumberOfColumns() const {
  return ncols;
}

inline
bool Value::isSymmetric() const {
  return symmetric;
}

inline
unsigned Value::getNumberOfGridDerivatives() const {
  return ngrid_der;
}

inline
double Value::getGridDerivative(const unsigned& n, const unsigned& j ) const {
  plumed_dbg_assert( hasDeriv && n*(1+ngrid_der) + 1 + j < data.size() );
  return data[n*(1+ngrid_der) + 1 + j];
}

inline
void Value::addGridDerivatives( const unsigned& n, const unsigned& j, const double& val ) {
  plumed_dbg_assert( hasDeriv && n*(1+ngrid_der) + 1 + j < data.size() );
  data[n*(1+ngrid_der) + 1 + j] += val;
}

inline
void Value::setGridDerivatives( const unsigned& n, const unsigned& j, const double& val ) {
  plumed_dbg_assert( hasDeriv && n*(1+ngrid_der) + 1 + j < data.size() );
  data[n*(1+ngrid_der) + 1 + j] = val;
}

}
#endif

