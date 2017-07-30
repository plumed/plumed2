/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2017 The plumed team
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
#ifndef __PLUMED_symfunc_SymmetryFunctionBase_h
#define __PLUMED_symfunc_SymmetryFunctionBase_h

#include <vector>
#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "adjmat/AdjacencyMatrixBase.h"

namespace PLMD {
namespace symfunc {

class SymmetryFunctionBase : 
public ActionWithValue,
public ActionWithArguments
{
private:
  unsigned nderivatives;
  std::vector<double> forcesToApply;
protected:
  void addValueWithDerivatives();
  void addComponentWithDerivatives( const std::string& name );
  void addToValue( const unsigned& ival, const double& val, MultiValue& myvals ) const ;
  void addWeightDerivative( const unsigned& ival, const double& dval, MultiValue& myvals ) const ;
  void addVectorDerivatives( const unsigned& ival, const Vector& dvec, MultiValue& myvals ) const ;
  void updateDerivativeIndices( MultiValue& myvals ) const ;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandMatrix( const bool& components, const std::string& lab, const std::vector<std::string>& words,
                                         const std::map<std::string,std::string>& keys,
                                         std::vector<std::vector<std::string> >& actions );
  static void expandFunctions( const std::string& labout, const std::string& argin, 
                               const std::vector<std::string>& words,
                               const std::map<std::string,std::string>& keys,
                               std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit SymmetryFunctionBase(const ActionOptions&); 
  unsigned getNumberOfDerivatives() const ;
  void calculate(){} 
  void buildCurrentTaskList( std::vector<unsigned>& tflags ) const ;
  virtual void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const = 0;
  virtual void performTask( const unsigned& current, MultiValue& myvals ) const ;
  void apply();
};

inline
unsigned SymmetryFunctionBase::getNumberOfDerivatives() const {
  return nderivatives;
}

inline
void SymmetryFunctionBase::addToValue( const unsigned& ival, const double& val, MultiValue& myvals ) const {
  myvals.addValue( getPntrToOutput(ival)->getPositionInStream(), val );
}

inline
void SymmetryFunctionBase::addWeightDerivative( const unsigned& ival, const double& dval, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  unsigned ostrn = getPntrToOutput(ival)->getPositionInStream();
  unsigned my_weight = getPntrToArgument(0)->getPositionInStream(); 
  for(unsigned k=0;k<myvals.getNumberActive(my_weight);++k){
      unsigned kind=myvals.getActiveIndex(my_weight,k); 
      myvals.addDerivative( ostrn, arg_deriv_starts[0] + kind, dval*myvals.getDerivative( my_weight, kind ) );
  }
}

inline
void SymmetryFunctionBase::addVectorDerivatives( const unsigned& ival, const Vector& dvec, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  unsigned ostrn = getPntrToOutput(ival)->getPositionInStream();
  unsigned my_x = getPntrToArgument(1)->getPositionInStream();
  for(unsigned k=0;k<myvals.getNumberActive(my_x);++k){
      unsigned kind=myvals.getActiveIndex(my_x,k);
      myvals.addDerivative( ostrn, arg_deriv_starts[1] + kind, dvec[0]*myvals.getDerivative( my_x, kind ) );
  }
  unsigned my_y = getPntrToArgument(2)->getPositionInStream();
  for(unsigned k=0;k<myvals.getNumberActive(my_y);++k){
      unsigned kind=myvals.getActiveIndex(my_y,k);
      myvals.addDerivative( ostrn, arg_deriv_starts[2] + kind, dvec[1]*myvals.getDerivative( my_y, kind ) );
  }
  unsigned my_z = getPntrToArgument(3)->getPositionInStream();
  for(unsigned k=0;k<myvals.getNumberActive(my_z);++k){
      unsigned kind=myvals.getActiveIndex(my_z,k);
      myvals.addDerivative( ostrn, arg_deriv_starts[3] + kind, dvec[2]*myvals.getDerivative( my_z, kind ) );
  }
}

}
}
#endif

