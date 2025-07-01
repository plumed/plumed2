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
#ifndef __PLUMED_function_FunctionOfScalar_h
#define __PLUMED_function_FunctionOfScalar_h

#include "Function.h"
#include "FunctionSetup.h"

namespace PLMD {
namespace function {

/**
\ingroup INHERIT
This is the abstract base class to use for implementing new CV function, within it there is
\ref AddingAFunction "information" as to how to go about implementing a new function.
*/

template <class T>
class FunctionOfScalar : public Function {
private:
/// Set equal to one if we are doing evaluateGridFunction
  unsigned argstart;
/// The function that is being computed
  T myfunc;
public:
  explicit FunctionOfScalar(const ActionOptions&);
  virtual ~FunctionOfScalar() {}
/// Get the label to write in the graph
  std::string writeInGraph() const override ;
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  void calculate() override;
  static void registerKeywords(Keywords&);
};

template <class T>
void FunctionOfScalar<T>::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  std::string name = keys.getDisplayName();
  std::size_t und=name.find("_SCALAR");
  keys.setDisplayName( name.substr(0,und) );
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  T::registerKeywords( keys );
  if( keys.getDisplayName()=="SUM" ) {
    keys.setValueDescription("scalar","the sum of all the input arguments");
  } else if( keys.getDisplayName()=="MEAN" ) {
    keys.setValueDescription("scalar","the mean of all the input arguments");
  } else if( keys.getDisplayName()=="EVALUATE_FUNCTION_FROM_GRID" ) {
    keys.addInputKeyword("compulsory","ARG","scalar/grid","");
  }
}

template <class T>
FunctionOfScalar<T>::FunctionOfScalar(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  argstart(0) {
  // Check if first argument is grid
  if( getPntrToArgument(0)->getRank()>0 && getPntrToArgument(0)->hasDerivatives() ) {
    argstart=1;
  }
  // Create the values to hold the output
  std::vector<std::size_t> shape;
  FunctionData<T>::setup( myfunc, keywords.getOutputComponents(), shape, true, this );
}

template <class T>
std::string FunctionOfScalar<T>::writeInGraph() const {
  std::size_t und = getName().find_last_of("_");
  return getName().substr(0,und);
}

template <class T>
std::string FunctionOfScalar<T>::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  if( getName().find("SORT")==std::string::npos ) {
    return ActionWithValue::getOutputComponentDescription( cname, keys );
  }
  return "the " + cname + "th largest of the input scalars";
}

template <class T>
void FunctionOfScalar<T>::calculate() {
  std::vector<double> args( getNumberOfArguments() - argstart );
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    args[i-argstart]=getPntrToArgument(i)->get();
  }
  std::vector<double> vals( getNumberOfComponents() ), deriv( getNumberOfComponents()*args.size() );
  auto funcout = FunctionOutput::create( getNumberOfComponents(),
                                         vals.data(),
                                         args.size(),
                                         deriv.data() );
  T::calc( myfunc,
           doNotCalculateDerivatives(),
           View<const double>(args.data(),args.size()),
           funcout );
  for(unsigned i=0; i<vals.size(); ++i) {
    copyOutput(i)->set(vals[i]);
  }
  if( doNotCalculateDerivatives() ) {
    return;
  }
  unsigned dstart = 0;
  for(unsigned i=0; i<argstart; ++i) {
    dstart += getPntrToArgument(i)->getNumberOfStoredValues();
  }

  for(unsigned i=0; i<vals.size(); ++i) {
    Value* val = getPntrToComponent(i);
    for(unsigned j=0; j<args.size(); ++j) {
      setDerivative( val, dstart+j, funcout.derivs[i][j] );
    }
  }
}

}
}
#endif
