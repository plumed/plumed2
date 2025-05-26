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
#ifndef __PLUMED_function_FunctionSetup_h
#define __PLUMED_function_FunctionSetup_h

#include "core/ActionWithValue.h"
#include "core/ActionWithArguments.h"
#include "core/ActionWithVector.h"

namespace PLMD {
namespace function {

class FunctionOptions {
public:
  /// Is the derivative zero if the value is zero
  bool derivativeZeroIfValueIsZero = false;
/// Are multiple components registered with a single name as in SphericalHarmonic and Moments
  std::vector<std::string> multipleValuesForEachRegisteredComponent;
};

template <class T>
class FunctionData {
public:
  /// Set equal to one if we are doing EvaluateGridFunction
  unsigned argstart = 0;
  // Number of scalars that appear in input 
  unsigned nscalars = 0;
  T f;
  FunctionData<T>& operator=( const FunctionData<T>& m ) {
    argstart = m.argstart;
    nscalars = m.nscalars;
    f = m.f;
    return *this;
  }
  // This is for setting up the functions
  static void setup( T& myfunc, const std::vector<std::string>& components, const std::vector<std::size_t>& shape, bool hasderiv, ActionWithValue* action );
};

template <class T>
void FunctionData<T>::setup( T& myfunc, const std::vector<std::string>& components, const std::vector<std::size_t>& shape, bool hasderiv, ActionWithValue* action ) {
  ActionWithArguments* aarg = dynamic_cast<ActionWithArguments*>( action );
  plumed_assert( aarg );
  FunctionOptions options;
  T::read( myfunc, aarg, options );

  if( action->keywords.getDisplayName()=="SORT" ) {
    for(unsigned j=0; j<aarg->getNumberOfArguments(); ++j) {
      std::string num;
      Tools::convert( j+1, num );
      if( hasderiv ) {
        action->addComponentWithDerivatives( num, shape );
      } else {
        action->addComponent( num, shape );
      }
    }
  } else {
    for(unsigned i=0; i<components.size(); ++i) {
      if( options.multipleValuesForEachRegisteredComponent.size()>0 ) {
        for(unsigned j=0; j<options.multipleValuesForEachRegisteredComponent.size(); ++j) {
          if( hasderiv ) {
            action->addComponentWithDerivatives( components[i] + options.multipleValuesForEachRegisteredComponent[j], shape );
          } else {
            action->addComponent( components[i] + options.multipleValuesForEachRegisteredComponent[j], shape );
          }
        }
      } else if( components[i]==".#!value" ) {
        if( hasderiv ) {
          action->addValueWithDerivatives( shape );
        } else {
          action->addValue( shape );
        }
      } else if( components[i].find_first_of("_")!=std::string::npos ) {
        ActionWithArguments* aarg = dynamic_cast<ActionWithArguments*>( action );
        plumed_assert( aarg );
        if( aarg->getNumberOfArguments()==1 && hasderiv ) {
          action->addValueWithDerivatives( shape );
        } else if( aarg->getNumberOfArguments()==1 ) {
          action->addValue( shape );
        } else {
          for(unsigned j=0; j<aarg->getNumberOfArguments(); ++j) {
            if( hasderiv ) {
              action->addComponentWithDerivatives( (aarg->getPntrToArgument(j))->getName() + components[i], shape );
            } else {
              action->addComponent( (aarg->getPntrToArgument(j))->getName() + components[i], shape );
            }
          }
        }
      } else {
        if( hasderiv ) {
          action->addComponentWithDerivatives( components[i], shape );
        } else {
          action->addComponent( components[i], shape );
        }
      }
    }
  }
  if( options.derivativeZeroIfValueIsZero ) {
    for(int i=0; i<action->getNumberOfComponents(); ++i) {
      (action->copyOutput(i))->setDerivativeIsZeroWhenValueIsZero();
    }
  }
  if( action->keywords.exists("PERIODIC") ) {
    plumed_assert( action->keywords.getDisplayName()!="DIFFERENCE" );
    std::vector<std::string> period;
    action->parseVector("PERIODIC",period);
    if( period.size()==1 ) {
      if( period[0]!="NO") {
        action->error("input to PERIODIC keyword does not make sense");
      }
      action->setNotPeriodic();
      return;
    } else if( period.size()!=2 ) {
      action->error("input to PERIODIC keyword does not make sense");
    }
    for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
      (action->copyOutput(i))->setDomain( period[0], period[1] );
    }
  } else if( action->keywords.getDisplayName()=="DIFFERENCE" ) {
    ActionWithArguments* aarg = dynamic_cast<ActionWithArguments*>( action );
    plumed_assert( aarg );
    if( aarg->getPntrToArgument(0)->isPeriodic() ) {
      std::string min, max;
      aarg->getPntrToArgument(0)->getDomain( min, max );
      action->setPeriodic(min,max);
    } else if( aarg->getPntrToArgument(1)->isPeriodic() ) {
      std::string min, max;
      aarg->getPntrToArgument(0)->getDomain( min, max );
      action->setPeriodic( min, max );
    } else {
      action->setNotPeriodic();
    }
  } else {
    for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
      (action->copyOutput(i))->setNotPeriodic();
    }
  }
}

class FunctionOutput {
public:
  unsigned nvals;
  View<double,helpers::dynamic_extent> values;
  unsigned nder;
  View2D<double,helpers::dynamic_extent,helpers::dynamic_extent> derivs;
  FunctionOutput( unsigned nv, double* v, unsigned na, double* d ):
    nvals(nv),
    values(v,nv),
    nder(na),
    derivs(d,nv,na) {
  }
};

}
}
#endif
