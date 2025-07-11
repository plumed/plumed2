/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "FunctionSetup.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION PIECEWISE
/*
Compute a piece wise straight line through its arguments that passes through a set of ordered control points.

This action can be used to calculate a piecewise linear function of an input argument such as the one given below:

$$
f(x) = \begin{cases}
10 & \textrm{if} \quad x<1 \\
10 + \frac{\pi - 10}{2-1}(x-1) & \textrm{if} \quad 1 \le x < 2 \\
\pi + \frac{10 - \pi}{3-2}(x-2) & \textrm{if} \quad 2 \le x \le 3 \\
10 & \textrm{otherwise}
\end{cases}
$$

The example shown below illustrates how one can use PLUMED to evaluate the function described above for the distance
between atom 1 and atom 2.

```plumed
dist1: DISTANCE ATOMS=1,10
pw: PIECEWISE POINT0=1,10 POINT1=2,PI POINT2=3,10 ARG=dist1
PRINT ARG=pw FILE=colvar
```

As you can see from the example above, the control points for the picewise function are passed using the `POINT0=...` `POINT1=...` syntax.
You can specify as many of these control points as you want.  These control points then serce as the $x_i$ and $y_i$ values in the following expression.

For variables less than the first
(greater than the last) point, the value of the first (last) point is used.

$$
\frac{y_{i+1}-y_i}{x_{i+1}-x_i}(s-x_i)+y_i ;  if x_i<s<x_{i+1}
$$

If the input value of $s$ is smaller than the lowest specified $x_i$ value then this action outputs the $y_i$ value that corresponds to the smallest of the input $x_i$ values.
Similarly, if the input value of $s$ is larger than the highest specified $x_i$ value the $y_i$ value that corresponds to the largest of the input $x_i$ values is output.

## Using multiple scalars in input

The following example illustrates what happens when multiple scalar arguments are passed to this action:

```plumed
dist1: DISTANCE ATOMS=1,10
dist2: DISTANCE ATOMS=2,11

ppww: PIECEWISE POINT0=1,10 POINT1=2,PI POINT2=3,10 ARG=dist1,dist2
PRINT ARG=ppww.dist1_pfunc,ppww.dist2_pfunc
```

In essence the piecewise function is applied to each of the input arguments in turn.  Hence, in the example above the PIECEWISE command outputs two values.  The first of
these `ppww.dist1_pfunc` is the result that is obtained when the piecewise function is applied to the argument `dist1`.  The second is then the result that is obtained when the
piecewise function is applied to the argument `dist2`.

## Non rank zero arguments

This argument currently cannot accept non-rank zero arguments.  However, it would be extremely straightforward to add functionality to ensure that if a PIECEWISE command receives
a vector, matrix or function on grid in input it will output a vector, matrix or function on grid that is obtained by applying the piecewise function elementwise to each of the
elements of the input vector, matrix or function.

*/
//+ENDPLUMEDOC

class Piecewise {
public:
  std::vector<std::pair<double,double> > points;
  static void registerKeywords(Keywords& keys);
  static void read( Piecewise& func, ActionWithArguments* action, FunctionOptions& options );
  static void calc( const Piecewise& func, bool noderiv, View<const double> args, FunctionOutput& funcout );
};


typedef FunctionShortcut<Piecewise> PiecewiseShortcut;
PLUMED_REGISTER_ACTION(PiecewiseShortcut,"PIECEWISE")
typedef FunctionOfScalar<Piecewise> ScalarPiecewise;
PLUMED_REGISTER_ACTION(ScalarPiecewise,"PIECEWISE_SCALAR")

void Piecewise::registerKeywords(Keywords& keys) {
  keys.add("numbered","POINT","This keyword is used to specify the various points in the function above.");
  keys.reset_style("POINT","compulsory");
  keys.addOutputComponent("_pfunc","default","scalar","one or multiple instances of this quantity can be referenced elsewhere "
                          "in the input file.  These quantities will be named with the arguments of the "
                          "function followed by the character string _pfunc.  These quantities tell the "
                          "user the values of the piece wise functions of each of the arguments.");
}

void Piecewise::read( Piecewise& func, ActionWithArguments* action, FunctionOptions& options ) {
  for(int i=0;; i++) {
    std::vector<double> pp;
    if(!action->parseNumberedVector("POINT",i,pp) ) {
      break;
    }
    if(pp.size()!=2) {
      action->error("points should be in x,y format");
    }
    func.points.push_back(std::pair<double,double>(pp[0],pp[1]));
    if(i>0 && func.points[i].first<=func.points[i-1].first) {
      action->error("points abscissas should be monotonously increasing");
    }
  }

  for(unsigned i=0; i<action->getNumberOfArguments(); i++) {
    if(action->getPntrToArgument(i)->isPeriodic()) {
      action->error("Cannot use PIECEWISE on periodic arguments");
    }
  }
  action->log.printf("  on points:");
  for(unsigned i=0; i<func.points.size(); i++) {
    action->log.printf("   (%f,%f)",func.points[i].first,func.points[i].second);
  }
  action->log.printf("\n");
}

void Piecewise::calc( const Piecewise& func, bool noderiv, const View<const double> args, FunctionOutput& funcout ) {
  for(unsigned i=0; i<args.size(); i++) {
    unsigned p=0;
    for(; p<func.points.size(); p++) {
      if(args[i]<func.points[p].first) {
        break;
      }
    }
    if(p==0) {
      funcout.values[i]=func.points[0].second;
      if( !noderiv ) {
        funcout.derivs[i][i]=0.0;
      }
    } else if(p==func.points.size()) {
      funcout.values[i]=func.points[func.points.size()-1].second;
      if( !noderiv ) {
        funcout.derivs[i][i]=0.0;
      }
    } else {
      double m=(func.points[p].second-func.points[p-1].second) / (func.points[p].first-func.points[p-1].first);
      funcout.values[i]=m*(args[i]-func.points[p-1].first)+func.points[p-1].second;
      if( !noderiv ) {
        funcout.derivs[i][i]=m;
      }
    }
  }
}

}
}


