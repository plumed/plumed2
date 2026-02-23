/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#ifndef __PLUMED_symfunc_Fccubic_h
#define __PLUMED_symfunc_Fccubic_h
#include "function/FunctionSetup.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace symfunc {

template<typename T=double>
struct Fccubic {
  using precision=T;
  precision alpha;
  precision a1;
  precision b1;
  static void registerKeywords( Keywords& keys );
  static void read( Fccubic& func,
                    ActionWithArguments* action,
                    function::FunctionOptions& funcout );
  static void calc( const Fccubic& func,
                    bool noderiv,
                    View<const precision> args,
                    function::FunctionOutput& funcout );

#ifdef __PLUMED_HAS_OPENACC
  void toACCDevice() const {
#pragma acc enter data copyin(this[0:1], alpha,a1,b1)
  }
  void removeFromACCDevice() const {
#pragma acc exit data delete(b1,a1,alpha,this[0:1])
  }
#endif // __PLUMED_HAS_OPENACC
};

template <typename T>
void Fccubic<T>::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","ALPHA","3.0",
           "The alpha parameter of the angular function");
  keys.setValueDescription("matrix",
                           "a function that measures the similarity with an fcc environment");
}

template <typename T>
void Fccubic<T>::read( Fccubic& func,
                       ActionWithArguments* action,
                       function::FunctionOptions& funcout ) {
  // Scaling factors such that '1' corresponds to fcc lattice
  // and '0' corresponds to isotropic (liquid)
  action->parse("ALPHA",func.alpha);
  func.a1 = 80080. / (2717. + 16*func.alpha);
  func.b1 = 16.*(func.alpha-143)/(2717+16*func.alpha);
  action->log.printf("  setting alpha paramter equal to %f \n",func.alpha);
}

template <typename T>
void Fccubic<T>::calc( const Fccubic& func,
                       bool noderiv,
                       const View<const precision> args,
                       function::FunctionOutput& funcout ) {
  precision x2 = args[0]*args[0];
  precision x4 = x2*x2;

  precision y2 = args[1]*args[1];
  precision y4 = y2*y2;

  precision z2 = args[2]*args[2];
  precision z4 = z2*z2;

  precision d2 = x2 + y2 + z2;
  if( d2 < epsilon ) {
    d2 = 1;
  }
  precision r8 = pow( d2, 4 );
  precision r12 = pow( d2, 6 );

  precision tmp = ((x4*y4)+(x4*z4)+(y4*z4))/r8-func.alpha*x4*y4*z4/r12;

  precision t0 = (x2*y4+x2*z4)/r8-func.alpha*x2*y4*z4/r12;
  precision t1 = (y2*x4+y2*z4)/r8-func.alpha*y2*x4*z4/r12;
  precision t2 = (z2*x4+z2*y4)/r8-func.alpha*z2*x4*y4/r12;
  precision t3 = (2*tmp-func.alpha*x4*y4*z4/r12)/d2;

  if( !noderiv ) {
    funcout.derivs[0][0]=4*func.a1*args[0]*(t0-t3);
    funcout.derivs[0][1]=4*func.a1*args[1]*(t1-t3);
    funcout.derivs[0][2]=4*func.a1*args[2]*(t2-t3);
  }

  // Set the value and the derivatives
  funcout.values[0] = (func.a1*tmp+func.b1);
}

}
}

#endif //__PLUMED_symfunc_Fccubic_h
