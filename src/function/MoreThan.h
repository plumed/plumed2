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
#ifndef __PLUMED_function_MoreThan_h
#define __PLUMED_function_MoreThan_h
#include "FunctionSetup.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionOfMatrix.h"
#include "core/ActionRegister.h"
#include "tools/SwitchingFunction.h"

#include <cmath>

namespace PLMD {
namespace function {

class MoreThan {
public:
  bool squared;
#ifdef __PLUMED_HAS_OPENACC
  SwitchingFunctionAccelerable switchingFunction;
#else
  SwitchingFunction switchingFunction;
#endif //__PLUMED_HAS_OPENACC
  static void registerKeywords( Keywords& keys );
  static void read( MoreThan& func,
                    ActionWithArguments* action,
                    FunctionOptions& options );
  static void calc( const MoreThan& func,
                    bool noderiv,
                    const View<const double>& args,
                    FunctionOutput& funcout );
#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1],squared)
    switchingFunction.toACCDevice();
  }
  void removeFromACCDevice() const {
    switchingFunction.removeFromACCDevice();
#pragma acc exit data delete(squared,this[0:1])
  }
#endif //__PLUMED_HAS_OPENACC
};

void MoreThan::registerKeywords(Keywords& keys) {
  keys.add("compulsory","NN","6","The n parameter of the switching function ");
  keys.add("compulsory","MM","0","The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","D_0","0.0","The d_0 parameter of the switching function");
  keys.add("compulsory","R_0","The r_0 parameter of the switching function");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.addFlag("SQUARED",false,"is the input quantity the square of the value that you would like to apply the switching function to");
  keys.setValueDescription("scalar/vector/matrix","a function that is one if the if the input is more than a threshold");
}

void MoreThan::read( MoreThan& func, ActionWithArguments* action, FunctionOptions& options ) {
  options.derivativeZeroIfValueIsZero = true;
  if( action->getNumberOfArguments()!=1 ) {
    ActionWithVector* av = dynamic_cast<ActionWithVector*>( action );
    if( !av || (av && action->getNumberOfArguments()-av->getNumberOfMasks()!=1) ) {
      action->error("should only be one argument to less_than actions");
    }
  }
  if( action->getPntrToArgument(0)->isPeriodic() ) {
    action->error("cannot use this function on periodic functions");
  }


  std::string errors;
  std::string sfinput;
  action->parse("SWITCH", sfinput);
  if(sfinput.length()>0) {
    func.switchingFunction.set(sfinput, errors);
    if( errors.length()!=0 ) {
      action->error("problem reading SWITCH keyword : " + errors );
    }
  } else {
    int nn=6;
    int mm=0;
    double d0=0.0;
    double r0=0.0;
    action->parse("R_0",r0);
    if(r0<=0.0) {
      action->error("R_0 should be explicitly specified and positive");
    }
    action->parse("D_0",d0);
    action->parse("NN",nn);
    action->parse("MM",mm);
    func.switchingFunction.set(nn,mm,r0,d0);
  }
  action->log<<"  using switching function with cutoff "<<func.switchingFunction.description()<<"\n";
  action->parseFlag("SQUARED",func.squared);
  if( func.squared ) {
    action->log<<"  input quantity is square of quantity that switching function acts upon\n";
  }
}

void MoreThan::calc( const MoreThan& func, bool noderiv, const View<const double>& args, FunctionOutput& funcout ) {
  // the presence of NDEBUG seems to be ignored by nvcc...
  // plumed_dbg_assert( args.size()==1 );
  double d;
  if( func.squared ) {
    funcout.values[0] = 1.0 - func.switchingFunction.calculateSqr( args[0], d );
  } else {
    funcout.values[0] = 1.0 - func.switchingFunction.calculate( args[0], d );
  }
  if( !noderiv ) {
    funcout.derivs[0][0] = -args[0]*d;
  }
}

}
}
#endif //__PLUMED_function_MoreThan_h
