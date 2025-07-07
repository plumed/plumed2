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
#ifndef __PLUMED_function_FunctionWithSingleArgument_h
#define __PLUMED_function_FunctionWithSingleArgument_h

#include "core/ActionWithArguments.h"
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "FunctionSetup.h"

namespace PLMD {
namespace function {

template <class T>
class FunctionWithSingleArgument :
  public ActionWithValue,
  public ActionWithArguments {
private:
  T f;
  bool ismatrix;
public:
  static void registerKeywords(Keywords& keys);
  explicit FunctionWithSingleArgument(const ActionOptions&);
  /// Get the label to write in the graph
  std::string writeInGraph() const override {
    std::size_t und = getName().find_last_of("_");
    return getName().substr(0,und);
  }
  unsigned getNumberOfDerivatives() override ;
  void calculate() override ;
  void apply() override ;
};

template <class T>
void FunctionWithSingleArgument<T>::registerKeywords(Keywords& keys) {
  Action::registerKeywords(keys);
  ActionWithValue::registerKeywords(keys);
  ActionWithArguments::registerKeywords(keys);
  std::string name = keys.getDisplayName();
  std::size_t und=name.find("_ONEARG");
  keys.setDisplayName( name.substr(0,und) );
  keys.addInputKeyword("compulsory","ARG","vector/matrix","the vector/matrix that is being used in input");
  T::registerKeywords( keys );
}

template <class T>
FunctionWithSingleArgument<T>::FunctionWithSingleArgument(const ActionOptions&ao):
  Action(ao),
  ActionWithValue(ao),
  ActionWithArguments(ao),
  ismatrix(false) {

  if( getNumberOfArguments()!=1 ) {
    error("should only be one argument in input for this action");
  }
  ismatrix = getPntrToArgument(0)->getRank()==2 && !getPntrToArgument(0)->hasDerivatives();

  std::vector<std::size_t> shape;
  FunctionData<T>::setup( f, keywords.getOutputComponents(), shape, true, this );
}

template <class T>
unsigned FunctionWithSingleArgument<T>::getNumberOfDerivatives() {
  return getPntrToArgument(0)->getNumberOfStoredValues();
}

template <class T>
void FunctionWithSingleArgument<T>::calculate() {
  const Value* myarg = getPntrToArgument(0);
  std::vector<double> args( myarg->getNumberOfStoredValues() );
  if( ismatrix ) {
    unsigned nvals = 0;
    for(unsigned i=0; i<myarg->getShape()[0]; ++i) {
      nvals += myarg->getRowLength(i);
    }
    args.resize( nvals );
    unsigned k=0;
    unsigned ncols = myarg->getNumberOfColumns();
    for(unsigned i=0; i<myarg->getShape()[0]; ++i) {
      for(unsigned j=0; j<myarg->getRowLength(i); ++j) {
        args[k] = myarg->get(i*ncols + j, false);
        ++k;
      }
    }
  } else {
    for(unsigned i=0; i<args.size(); ++i) {
      args[i] = myarg->get( i );
    }
  }
  std::vector<double> vals( getNumberOfComponents() ), deriv( getNumberOfComponents()*args.size() );
  auto funcout = FunctionOutput::create( getNumberOfComponents(),
                                         vals.data(),
                                         args.size(),
                                         deriv.data() );
  T::calc( f,
           doNotCalculateDerivatives(),
           View<const double>(args.data(), args.size()),
           funcout );
  for(unsigned i=0; i<getNumberOfComponents(); ++i) {
    Value* myval=getPntrToComponent(i);
    myval->set( vals[i] );

    if( doNotCalculateDerivatives() ) {
      continue ;
    }

    for(unsigned j=0; j<args.size(); ++j) {
      myval->setDerivative( j, funcout.derivs[i][j] );
    }
  }
}

template <class T>
void FunctionWithSingleArgument<T>::apply() {
  if( !getPntrToComponent(0)->forcesWereAdded() ) {
    return ;
  }

  Value* myarg = getPntrToArgument(0);
  if( ismatrix ) {
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      unsigned k=0;
      const Value* myval = getPntrToComponent(i);
      double force = myval->getForce();
      unsigned ncols = myarg->getNumberOfColumns();
      for(unsigned n=0; n<myarg->getShape()[0]; ++n) {
        for(unsigned j=0; j<myarg->getRowLength(i); ++j) {
          myarg->addForce( n*ncols+j, force*myval->getDerivative(k), false );
          ++k;
        }
      }
    }
  } else {
    for(unsigned i=0; i<getNumberOfComponents(); ++i) {
      const Value* myval = getPntrToComponent(i);
      double force = myval->getForce();
      for(unsigned j=0; j<myarg->getNumberOfStoredValues(); ++j) {
        myarg->addForce( j, force*myval->getDerivative(j) );
      }
    }
  }
}

}
}
#endif
