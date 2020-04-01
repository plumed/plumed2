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

class ActionShortcut;

namespace symfunc {

class SymmetryFunctionBase :
  public ActionWithValue,
  public ActionWithArguments
{
private:
  bool usecols;
  unsigned nderivatives;
  std::vector<double> forcesToApply;
protected:
  bool done_with_matrix_comput;
  void addValueWithDerivatives();
  void addComponentWithDerivatives( const std::string& name );
  void addToValue( const unsigned& ival, const double& val, MultiValue& myvals ) const ;
  void addWeightDerivative( const unsigned& ival, const double& dval, MultiValue& myvals ) const ;
  void addVectorDerivatives( const unsigned& ival, const Vector& dvec, MultiValue& myvals ) const ;
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandMatrix( const bool& components, const std::string& lab, const std::string& sp_str,
                            const std::string& spa_str, const std::string& spb_str, ActionShortcut* action );
  static void createSymmetryFunctionObject( const std::string& lab, const std::string& name, const bool& iscoord, const bool& norm, ActionShortcut* action );
  static void registerKeywords( Keywords& keys );
  explicit SymmetryFunctionBase(const ActionOptions&);
  bool mustBeTreatedAsDistinctArguments() const ;
  virtual unsigned getNumberOfDerivatives() const ;
  virtual void calculate(){}
  virtual void compute( const double& weight, const Vector& vec, MultiValue& myvals ) const = 0;
  virtual void updateDerivativeIndices( MultiValue& myvals ) const ;
  virtual void performTask( const unsigned& current, MultiValue& myvals ) const ;
  virtual void computeSymmetryFunction( const unsigned& current, MultiValue& myvals ) const ;
  void apply();
};

inline
bool SymmetryFunctionBase::mustBeTreatedAsDistinctArguments() const {
  return true;
}

inline
unsigned SymmetryFunctionBase::getNumberOfDerivatives() const {
  if( usecols ) return 0;
  return nderivatives;
}

inline
void SymmetryFunctionBase::addToValue( const unsigned& ival, const double& val, MultiValue& myvals ) const {
  if( usecols ) {
    unsigned col_stash_index = myvals.getSecondTaskIndex(); if( col_stash_index>=getFullNumberOfTasks() ) col_stash_index -= getFullNumberOfTasks();
    myvals.stashMatrixElement( getPntrToOutput(ival)->getPositionInMatrixStash(), col_stash_index, val );
  } else myvals.addValue( getPntrToOutput(ival)->getPositionInStream(), val );
}

inline
void SymmetryFunctionBase::addWeightDerivative( const unsigned& ival, const double& dval, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  if( done_with_matrix_comput ) {
    unsigned ostrn = getPntrToOutput(ival)->getPositionInStream();
    unsigned my_weight = getPntrToArgument(0)->getPositionInStream();
    for(unsigned k=0; k<myvals.getNumberActive(my_weight); ++k) {
      unsigned kind=myvals.getActiveIndex(my_weight,k);
      myvals.addDerivative( ostrn, arg_deriv_starts[0] + kind, dval*myvals.getDerivative( my_weight, kind ) );
    }
  } else {
    unsigned index = ival*getPntrToArgument(0)->getShape()[1]+myvals.getSymfuncTemporyIndex();
    unsigned my_w = getPntrToArgument(0)->getPositionInStream();
    myvals.getSymfuncTemporyDerivatives(my_w)[index]+=dval;
  }
}

inline
void SymmetryFunctionBase::addVectorDerivatives( const unsigned& ival, const Vector& dvec, MultiValue& myvals ) const {
  if( doNotCalculateDerivatives() ) return;
  if( done_with_matrix_comput ) {
    unsigned ostrn = getPntrToOutput(ival)->getPositionInStream();
    for(unsigned i=0;i<getNumberOfArguments()-1;++i) {
        unsigned my_x = getPntrToArgument(1+i)->getPositionInStream();
        for(unsigned k=0; k<myvals.getNumberActive(my_x); ++k) {
          unsigned kind=myvals.getActiveIndex(my_x,k);
          myvals.addDerivative( ostrn, arg_deriv_starts[1+i] + kind, dvec[i]*myvals.getDerivative( my_x, kind ) );
        }
    }
  } else {
    unsigned index = ival*getPntrToArgument(0)->getShape()[1]+myvals.getSymfuncTemporyIndex();
    for(unsigned i=0;i<getNumberOfArguments()-1;++i) {
        unsigned my_x = getPntrToArgument(1+i)->getPositionInStream();
        myvals.getSymfuncTemporyDerivatives(my_x)[index]+=dvec[i];
    }
  }
}

}
}
#endif

