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
#include "tools/Matrix.h"

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
/// The function that is being computed
  T myfunc;
/// Are we on the first step
  bool firststep;
public:
  explicit FunctionOfScalar(const ActionOptions&);
  virtual ~FunctionOfScalar() {}
/// Get the label to write in the graph
  std::string writeInGraph() const override {
    return myfunc.getGraphInfo( getName() );
  }
  std::string getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const override ;
  void calculate() override;
  static void registerKeywords(Keywords&);
  void turnOnDerivatives() override;
};

template <class T>
void FunctionOfScalar<T>::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  std::string name = keys.getDisplayName();
  std::size_t und=name.find("_SCALAR");
  keys.setDisplayName( name.substr(0,und) );
  keys.add("hidden","NO_ACTION_LOG","suppresses printing from action on the log");
  T tfunc;
  tfunc.registerKeywords( keys );
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
  firststep(true) {
  myfunc.read( this );
  // Get the names of the components
  std::vector<std::string> components( keywords.getOutputComponents() );
  // Create the values to hold the output
  std::vector<std::string> str_ind( myfunc.getComponentsPerLabel() );
  for(unsigned i=0; i<components.size(); ++i) {
    if( str_ind.size()>0 ) {
      std::string compstr = components[i];
      if( compstr==".#!value" ) {
        compstr = "";
      }
      for(unsigned j=0; j<str_ind.size(); ++j) {
        addComponentWithDerivatives( compstr + str_ind[j] );
      }
    } else if( components[i]==".#!value" ) {
      addValueWithDerivatives();
    } else if( components[i].find_first_of("_")!=std::string::npos ) {
      if( getNumberOfArguments()==1 ) {
        addValueWithDerivatives();
      } else {
        for(unsigned j=0; j<getNumberOfArguments(); ++j) {
          addComponentWithDerivatives( getPntrToArgument(j)->getName() + components[i] );
        }
      }
    } else {
      addComponentWithDerivatives( components[i] );
    }
  }
  // Set the periodicities of the output components
  myfunc.setPeriodicityForOutputs( this );
  myfunc.setPrefactor( this, 1.0 );
}

template <class T>
std::string FunctionOfScalar<T>::getOutputComponentDescription( const std::string& cname, const Keywords& keys ) const {
  if( getName().find("SORT")==std::string::npos ) {
    return ActionWithValue::getOutputComponentDescription( cname, keys );
  }
  return "the " + cname + "th largest of the input scalars";
}

template <class T>
void FunctionOfScalar<T>::turnOnDerivatives() {
  if( !myfunc.derivativesImplemented() ) {
    error("derivatives have not been implemended for " + getName() );
  }
  ActionWithValue::turnOnDerivatives();
}

template <class T>
void FunctionOfScalar<T>::calculate() {
  if( firststep ) {
    myfunc.setup( this );
    firststep=false;
  }
  unsigned argstart = myfunc.getArgStart();
  std::vector<double> args( getNumberOfArguments() - argstart );
  for(unsigned i=argstart; i<getNumberOfArguments(); ++i) {
    args[i-argstart]=getPntrToArgument(i)->get();
  }
  std::vector<double> vals( getNumberOfComponents() );
  Matrix<double> derivatives( getNumberOfComponents(), args.size() );
  myfunc.calc( this, args, vals, derivatives );
  for(unsigned i=0; i<vals.size(); ++i) {
    copyOutput(i)->set(vals[i]);
  }
  if( doNotCalculateDerivatives() ) {
    return;
  }

  for(unsigned i=0; i<vals.size(); ++i) {
    Value* val = getPntrToComponent(i);
    for(unsigned j=0; j<args.size(); ++j) {
      setDerivative( val, j, derivatives(i,j) );
    }
  }
}

}
}
#endif
