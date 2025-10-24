/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2024 The plumed team
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

#include "CLTool.h"
#include "core/CLToolRegister.h"
#include "tools/Tools.h"
#include "tools/SwitchingFunction.h"
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>

//+PLUMEDOC TOOLS plotswitch
/*
plotswitch is a tool that takes a the input of a switching function and tabulates the output on the terminal

The tabulated data is compatible with gnuplot and numpy.loadtxt so you can use this tool to plot any functions that you plan to use with commands like
[LESS_THAN](LESS_THAN.md) or [MORE_THAN](MORE_THAN.md).

Without options plotswitch will tabulate 50 points between 0 and `R_0`, and then continue in tabulating points with the same step until the value of `D_MAX` is reached. If
`D_MAX is unset then a value of `2R_0` is used in place of `D_MAX`.

Without options plotswitch will tabulate data calling calculateSqr, since should be the most used option within the various colvars

Note that if `R_0` happens to be between "from" and "to" the number of steps may not be exacly the number requested as this command forces the value of the function at r0 to be calculated.

The various `--rational**` options use the special set option for the rational, like in COORDINATION.

## Examples

If this command is called without options like this:

```plumed
plumed plotswitch > plot.dat
```

The function calculated is:

$$
f(x) = \frac{1 - x^6}{1 - x^{12}}
$$

If you want to use a function that is differnt to this one you use the `--switch` keyword as shown below:

```plumed
plumed plotswitch --switch="RATIONAL NN=5 MM=9 R_0=1.3" --from=1.29999 --to=1.30001 --steps=100 > plot.dat
```

The `--switch` keyword here takes the input for switching functions that is discussed in the documentation for
[LESS_THAN](LESS_THAN.md).  Notice also that if you use this command with an older plumed version you will see a discontinuity in dfunc at around 1.3
(if you use gnuplot with "p 'plot.dat' u 1:3 w l t 'dfunc', 'plot.dat' u 1:2 w l axis x1y2 t 'res'")

The following example shows another way of generating the `plot.dat` that is output by the command above:

```plumed
plumed plotswitch --rationalR_0=1.3 --rationalNN=5 --rationalMM=9 --rationalD_0=0 --from=1.29999 --to=1.30001 --steps=100 > plot.dat
```

As with [LESS_THAN](LESS_THAN.md), there is a alternative to `--switch` that can be used for sepcifying the parameters of RATIONAL switching function.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace cltools {
class SwitchingPlotter : public CLTool {
public:
  explicit SwitchingPlotter(const CLToolOptions&);
  static void registerKeywords( Keywords&  );

  int main( FILE*, FILE*, Communicator& ) override;

};
PLUMED_REGISTER_CLTOOL(SwitchingPlotter,"plotswitch")

void SwitchingPlotter::registerKeywords( Keywords& keys ) {
  CLTool::registerKeywords( keys );
  keys.add("compulsory","--switch",
           "RATIONAL NN=6 R_0=1.0","the input to pass to the switching function,"
           " please remeber the quotes");
  keys.add("compulsory","--steps","50",
           "the number of steps between 0 and R_O, or in the specified interval");
  keys.add("compulsory","--from","-1",
           "the start of the interval, if negative will be set to 0");
  keys.add("compulsory","--to","-1",
           "the end of the interval, will be D_MAX or 2*R_0 if D_MAX is not set");
  keys.add("compulsory","--plotprecision","8",
           "the precision to use for the tabulated results");
  keys.add("compulsory","--rationalR_0","-1",
           "The r_0 parameter of the switching function, this will activate the "
           "--rational options, note that this will ignore the --switch option silently");
  keys.add("compulsory","--rationalNN","6",
           "The n parameter of the switching function");
  keys.add("compulsory","--rationalMM","0",
           "The m parameter of the switching function; 0 implies 2*NN");
  keys.add("compulsory","--rationalD_0","0.0",
           "The d_0 parameter of the switching function");
  keys.addFlag("--nosquare",false,"use calculate instead of calculateSqr");
  keys.add("compulsory","--centerrange","-1",
           "centers the visualization in R_0 in a range given epsilons times r_0"
           ", note that specifying this will overide all the other range options");
}

SwitchingPlotter::SwitchingPlotter(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=inputType::commandline;
}
int SwitchingPlotter::main( FILE*, FILE*, Communicator& ) {
  //collecting options:
  std::string swInput;
  parse("--switch",swInput);
  bool dontOptimize;
  parseFlag("--nosquare",dontOptimize);
  int Nsteps;
  parse("--steps",Nsteps);
  double lowerLimit;
  parse("--from",lowerLimit);
  double upperLimit;
  parse("--to",upperLimit);
  unsigned plotPrecision;
  parse("--plotprecision",plotPrecision);
  int rationalNN;
  parse("--rationalNN",rationalNN);
  int rationalMM;
  parse("--rationalMM",rationalMM);
  double rationalD_0;
  parse("--rationalD_0",rationalD_0);
  double rationalR_0;
  parse("--rationalR_0",rationalR_0);
  //this works only because we use lepton to parse the numbers
  double centerrange;
  parse("--centerrange",centerrange);
  //setting up the switching function
  PLMD::SwitchingFunction switchingFunction;
  if (rationalR_0>0) {
    switchingFunction.set(rationalNN,rationalMM,rationalR_0,rationalD_0);
  } else {
    std::string errors;
    switchingFunction.set(swInput,errors);
    if( errors.length()!=0 ) {
      error("problem reading SWITCH keyword : " + errors );
    }
  }

  //setting up the limits:
  const double r0 = switchingFunction.get_r0();
  const double dmax = switchingFunction.get_dmax();

  if (lowerLimit <0) {
    lowerLimit=0.0;
  }
  if (upperLimit < 0) {
    upperLimit = dmax;
    if (! (upperLimit < std::numeric_limits<double>::max())) {
      upperLimit = 2*r0;
    }
  }
  if(centerrange>0) {
    upperLimit=(1.0+centerrange*PLMD::epsilon)*r0;
    lowerLimit=(1.0-centerrange*PLMD::epsilon)*r0;
  }
  const double step = [=]() {
    if(r0 > lowerLimit && r0< upperLimit) {
      //this will make the step pass trough r0
      double interval = (r0-lowerLimit)/(upperLimit-lowerLimit);
      return  (r0-lowerLimit)/(interval *Nsteps);
    }
    return (upperLimit-lowerLimit)/double(Nsteps);
  }
  ();
  if (step <0.0) {
    error("I calculated a negative step");
  }

  //finally doing the job
  //descriptions starts with the values of "r_0"
  std::cout <<"#r val dfunc ( r_0="<<switchingFunction.description()<<")\n";
  double x=lowerLimit;
  while(x < upperLimit) {
    double dfunc=0.0;
    double res;
    if(dontOptimize) {
      res=switchingFunction.calculate(x,dfunc);
    } else {
      res=switchingFunction.calculateSqr(x*x,dfunc);
    }
    std::cout << std::setprecision(plotPrecision) << x << "\t"
              << std::setprecision(plotPrecision) << res << "\t"
              << std::setprecision(plotPrecision) << dfunc << '\n';
    x+=step;
  }
  return 0;
}

} //namespace cltools
} // namespace PLMD
