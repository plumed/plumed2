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

struct FunctionOptions {
  /// Is the derivative zero if the value is zero
  bool derivativeZeroIfValueIsZero = false;
/// Are multiple components registered with a single name as in SphericalHarmonic and Moments
  std::vector<std::string> multipleValuesForEachRegisteredComponent;
};

template <class T>
struct FunctionData {
  /// Set equal to one if we are doing EvaluateGridFunction
  unsigned argstart = 0;
  // Number of scalars that appear in input
  unsigned nscalars = 0;
  T f;
  // This is for setting up the functions
  static void setup( T& myfunc,
                     const std::vector<std::string>& components,
                     const std::vector<std::size_t>& shape,
                     bool hasderiv,
                     ActionWithValue* action );
#ifdef __PLUMED_USE_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1], argstart, nscalars)
    f.toACCDevice();
  }
  void removeFromACCDevice() const {
    f.removeFromACCDevice();
#pragma acc exit data delete(nscalars, argstart, this[0:1])
  }
#endif //__PLUMED_USE_OPENACC
};

template <class T>
void FunctionData<T>::setup( T& myfunc,
                             const std::vector<std::string>& components,
                             const std::vector<std::size_t>& shape,
                             bool hasderiv,
                             ActionWithValue* action ) {
  ActionWithArguments* aarg = action->castToActionWithArguments();
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
    for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
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

struct FunctionOutput {
  View<double> values;
  View2D<double> derivs;
  static FunctionOutput create( unsigned nvals,
                                double* vals,
                                unsigned ndev,
                                double* devs ) {
    return FunctionOutput{
      View<double>(vals,nvals),
      View2D<double>(devs,nvals,ndev)
    };
  }
};

} // namespace function
} // namespace PLMD
#endif
