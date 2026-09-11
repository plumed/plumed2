/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#ifndef __PLUMED_function_Between_h
#define __PLUMED_function_Between_h
#include "FunctionSetup.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "tools/HistogramBead.h"

#include <cmath>

namespace PLMD {
namespace function {

struct Between {
  HistogramBead hist{HistogramBead::KernelType::gaussian,0.0,1.0,0.5};
  static void registerKeywords( Keywords& keys );
  static void read( Between& func, ActionWithArguments* action,
                    FunctionOptions& options );
  static void calc( const Between& func,
                    bool noderiv,
                    const View<const double> args,
                    FunctionOutput& funcout );
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1])
    hist.toACCDevice();
  }
  void removeFromACCDevice() const {
    hist.removeFromACCDevice();
#pragma acc exit data delete(order, this[0:1])
  }
#endif // __PLUMED_HAS_OPENACC
};

void Between::registerKeywords(Keywords& keys) {
  keys.add("compulsory","LOWER","the lower boundary for this particular bin");
  keys.add("compulsory","UPPER","the upper boundary for this particular bin");
  keys.add("compulsory","SMEAR","0.5","the ammount to smear the Gaussian for each value in the distribution");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous function defined above. "
           "The following provides information on the \\ref histogrambead that are available. "
           "When this keyword is present you no longer need the LOWER, UPPER, SMEAR and KERNEL keywords.");
  keys.setValueDescription("scalar/vector/matrix","a function that is one if the input falls within a particular range and zero otherwise");
}

void Between::read( Between& func, ActionWithArguments* action, FunctionOptions& options ) {
  options.derivativeZeroIfValueIsZero = true;
  if( action->getNumberOfArguments()!=1 ) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( action );
    if( !av || (av && action->getNumberOfArguments()-av->getNumberOfMasks()!=1) ) {
      action->error("should only be one argument to less_than actions");
    }
  }

  std::string str_min, str_max;
  bool isPeriodic = action->getPntrToArgument(0)->isPeriodic();
  if( isPeriodic ) {
    action->getPntrToArgument(0)->getDomain( str_min, str_max );
  }
  std::string hinput;
  action->parse("SWITCH",hinput);
  if(hinput.length()==0) {
    std::string low, up, sme;
    action->parse("LOWER",low);
    action->parse("UPPER",up);
    action->parse("SMEAR",sme);
    hinput = "GAUSSIAN LOWER=" + low + " UPPER=" + up + " SMEAR=" + sme;
  }
  std::string errors;
  func.hist.set( hinput, errors );
  if( errors.size()!=0 ) {
    action->error( errors );
  }
  action->log.printf("  %s \n", func.hist.description().c_str() );

  if( !isPeriodic ) {
    func.hist.isNotPeriodic();
  } else {
    double min;
    double max;
    Tools::convert( str_min, min );
    Tools::convert( str_max, max );
    func.hist.isPeriodic( min, max );
  }
}

void Between::calc( const Between& func,
                    bool noderiv,
                    const View<const double> args,
                    FunctionOutput& funcout ) {
  // the presence of NDEBUG seems to be ignored by nvcc...
  // plumed_dbg_assert( args.size()==1 );
  double deriv;
  funcout.values[0] = func.hist.calculate( args[0], deriv );
  if( !noderiv ) {
    funcout.derivs[0][0] = deriv;
  }
}

}
}
#endif //__PLUMED_function_Between_h
