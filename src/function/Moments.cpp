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
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"
#include "FunctionOfVector.h"
#include "FunctionWithSingleArgument.h"
#include "core/ActionRegister.h"
#include "FunctionSetup.h"

#include <cmath>

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION MOMENTS
/*
Calculate central moments from the distribution of input quantities

This action takes a set of $N$ input arguments, $s_i$, and evaluates the $k$th central moment of the distribution of input arguments using:

$$
\mu_k = \frac{1}{N} \sum_{i=1}^N ( s_i - \langle s \rangle )^k \qquad \textrm{where} \qquad \langle s \rangle = \frac{1}{N} \sum_{i=1}^N s_i
$$

A single moments action can evaluate more than one central moment at once so, for example, the input below can be used to calculate the second
and third central moment for the distribution of the four input distances.

```plumed
d12:  DISTANCE ATOMS=1,2
d13:  DISTANCE ATOMS=1,3
d14:  DISTANCE ATOMS=1,4
d15:  DISTANCE ATOMS=1,5
mv: MOMENTS ARG=d12,d13,d14,d15 POWERS=2,3
PRINT ARG=mv.moment-2,mv.moment-3 FILE=colvar
```

Notice that you can also achieve the same result using the following input:

```plumed
d: DISTANCE ATOMS1=1,2 ATOMS2=1,3 ATOMS3=1,4 ATOMS4=1,5
mv: MOMENTS ARG=d POWERS=2,3
PRINT ARG=mv.moment-2,mv.moment-3 FILE=colvar
```

In this second case the four distances are passed to the MOMENTS action as a vector.  The MOMENTS action then outputs 2 components - the
two central moments that were requested.

These examples are representative the only two ways you can use this action.  In input it can accept either a list of scalars or a single vector.
It does not accept matrices or a list of vectors in input.

*/
//+ENDPLUMEDOC

class Moments {
public:
  bool isperiodic;
  double min, max, pfactor;
  double max_minus_min;
  double inv_max_minus_min;
  std::vector<int> powers;
  static void registerKeywords(Keywords& keys);
  static void read( Moments& func, ActionWithArguments* action, FunctionOptions& options );
  static void calc( const Moments& func, bool noderiv, View<const double> args, FunctionOutput& funcout );
};

typedef FunctionShortcut<Moments> MomentsShortcut;
PLUMED_REGISTER_ACTION(MomentsShortcut,"MOMENTS")
typedef FunctionWithSingleArgument<Moments> SingleargMoments;
PLUMED_REGISTER_ACTION(SingleargMoments,"MOMENTS_ONEARG")
typedef FunctionOfScalar<Moments> ScalarMoments;
PLUMED_REGISTER_ACTION(ScalarMoments,"MOMENTS_SCALAR")
typedef FunctionOfVector<Moments> VectorMoments;
PLUMED_REGISTER_ACTION(VectorMoments,"MOMENTS_VECTOR")

void Moments::registerKeywords(Keywords& keys) {
  keys.add("compulsory","POWERS","calculate the central moments of the distribution of collective variables. "
           "The \\f$m\\f$th central moment of a distribution is calculated using \\f$\\frac{1}{N} \\sum_{i=1}^N ( s_i - \\overline{s} )^m \\f$, where \\f$\\overline{s}\\f$ is "
           "the average for the distribution. The POWERS keyword takes a lists of integers as input or a range. Each integer is a value of \\f$m\\f$. The final "
           "calculated values can be referenced using moment-\\f$m\\f$.");
  keys.addOutputComponent("moment","default","scalar","the central moments of the distribution of values. The second central moment "
                          "would be referenced elsewhere in the input file using "
                          "<em>label</em>.moment-2, the third as <em>label</em>.moment-3, etc.");
}

void Moments::read( Moments& func, ActionWithArguments* action, FunctionOptions& options ) {
  func.isperiodic = (action->getPntrToArgument(0))->isPeriodic();
  if( func.isperiodic ) {
    std::string str_min, str_max;
    (action->getPntrToArgument(0))->getDomain( str_min, str_max );
    for(unsigned i=1; i<action->getNumberOfArguments(); ++i) {
      if( !(action->getPntrToArgument(i))->isPeriodic() ) {
        action->error("cannot mix periodic and non periodic variables when calculating moments");
      }
      std::string str_min2, str_max2;
      (action->getPntrToArgument(i))->getDomain( str_min2, str_max2);
      if( str_min!=str_min2 || str_max!=str_max2 ) {
        action->error("all input arguments should have same domain when calculating moments");
      }
    }
    Tools::convert(str_min,func.min);
    Tools::convert(str_max,func.max);
    func.max_minus_min = func.max - func.min;
    func.inv_max_minus_min = 1 / func.max_minus_min;
    func.pfactor = 2*pi*func.inv_max_minus_min;
  } else {
    for(unsigned i=1; i<action->getNumberOfArguments(); ++i) {
      if( (action->getPntrToArgument(i))->isPeriodic() ) {
        action->error("cannot mix periodic and non periodic variables when calculating moments");
      }
    }
  }

  action->parseVector("POWERS",func.powers);
  for(unsigned i=0; i<func.powers.size(); ++i) {
    if( func.powers[i]<2 ) {
      action->error("first central moment is zero do you really need to calculate that");
    }
    action->log.printf("  computing %dth central moment of distribution of input cvs \n", func.powers[i]);
    std::string num;
    Tools::convert(func.powers[i],num);
    options.multipleValuesForEachRegisteredComponent.push_back( "-" + num );
  }
}

void Moments::calc( const Moments& func, bool noderiv, const View<const double> args, FunctionOutput& funcout ) {
  double mean=0;
  double inorm = 1.0 / static_cast<double>( args.size() );
  if( func.isperiodic ) {
    double sinsum=0, cossum=0, val;
    for(unsigned i=0; i<args.size(); ++i) {
      val=func.pfactor*( args[i] - func.min );
      sinsum+=sin(val);
      cossum+=cos(val);
    }
    mean = 0.5 + atan2( inorm*sinsum, inorm*cossum ) / (2*pi);
    mean = func.min + (func.max-func.min)*mean;
  } else {
    for(unsigned i=0; i<args.size(); ++i) {
      mean+=args[i];
    }
    mean *= inorm;
  }

  if( func.isperiodic ) {
    for(unsigned npow=0; npow<func.powers.size(); ++npow) {
      double dev1=0;
      for(unsigned i=0; i<args.size(); ++i) {
        double s = func.inv_max_minus_min*( args[i] - mean );
        s = Tools::pbc(s);
        dev1+=pow( s*func.max_minus_min, func.powers[npow] - 1 );
      }
      dev1*=inorm;
      funcout.values[npow] = 0;
      double prefactor = func.powers[npow]*inorm;
      for(unsigned i=0; i<args.size(); ++i) {
        double tmp = func.inv_max_minus_min*( args[i] - mean );
        tmp = func.max_minus_min*Tools::pbc(tmp);
        funcout.values[npow] += inorm*pow( tmp, func.powers[npow] );
        if( !noderiv ) {
          funcout.derivs[npow][i] = prefactor*(pow( tmp, func.powers[npow] - 1 ) - dev1);
        }
      }
    }
  } else {
    for(unsigned npow=0; npow<func.powers.size(); ++npow) {
      double dev1=0;
      for(unsigned i=0; i<args.size(); ++i) {
        dev1+=pow( args[i] - mean, func.powers[npow] - 1 );
      }
      dev1*=inorm;
      funcout.values[npow] = 0;
      double prefactor = func.powers[npow]*inorm;
      for(unsigned i=0; i<args.size(); ++i) {
        double tmp=args[i] - mean;
        funcout.values[npow] += inorm*pow( tmp, func.powers[npow] );
        if( !noderiv ) {
          funcout.derivs[npow][i] = prefactor*(pow( tmp, func.powers[npow] - 1 ) - dev1);
        }
      }
    }
  }
}

}
}


