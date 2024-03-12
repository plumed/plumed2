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
#include <cstdio>
#include <string>
#include <iostream>
#include <limits>

//+PLUMEDOC TOOLS plotswitch
/*
plotswitch is a tool that takes a the input of a switching function and tabulates the output on the terminal

The tabulated data is compatible with gnuplot and numpy.loadtxt

Without options plotswitch will tabulat 50 points between 0 and R_0, and then continue in tabulating points with the same step until 2*R_0 or if D_MAX is set, D_MAX

Without options plotswitch will tabulate data calling calculateSqr, since should be the most used option within the various colvars

The various --rational** options use the special set option for the rational, like in COORDINATION.

\par Examples

Without option will plot the NN=6 MM=12 rational
\verbatim
plumed plotswitch > plot.dat
\endverbatim

\verbatim
plumed plotswitch --switch="RATIONAL NN=5 MM=9 R_0=1.3" --from=1.29999 --to=1.30001 --steps=100> plot.dat
\endverbatim
If you use this with a older plumed version you will see the discontinuity in dfunc around 1.3
(i use gnuplot with "p 'plot.dat' u 1:3 w l t 'dfunc', 'plot.dat' u 1:2 w l axis x1y2 t 'res'")
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
}

SwitchingPlotter::SwitchingPlotter(const CLToolOptions& co ):
  CLTool(co) {
  inputdata=commandline;
}
int SwitchingPlotter::main( FILE*, FILE*, Communicator& ) {
  std::string swInput;
  parse("--switch",swInput);
  bool dontOptimize;
  parseFlag("--nosquare",dontOptimize);
  int Nsteps;
  parse("--steps",Nsteps);
  double lowerLimit;
  parse("--from",lowerLimit);
  if (lowerLimit <0) {
    lowerLimit=0.0;
  }
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

  std::string errors;
  PLMD::SwitchingFunction switchingFunction;
  if (rationalR_0>0) {
    switchingFunction.set(rationalNN,rationalMM,rationalR_0,rationalD_0);
  } else {
    switchingFunction.set(swInput,errors);

    if( errors.length()!=0 ) {
      error("problem reading SWITCH keyword : " + errors );
    }
  }
  // const double d0 = switchingFunction.get_d0();
  const double r0 = switchingFunction.get_r0();
  const double dmax = switchingFunction.get_dmax();
  const double step = [&]() {
    if (upperLimit < 0) {
      upperLimit = dmax;
      if (! (upperLimit < std::numeric_limits<double>::max())) {
        upperLimit = 2*r0;
        return r0/double(Nsteps);
      }
    }
    return (upperLimit-lowerLimit)/double(Nsteps);
  }();
  //descriptions starts with the values of "r_0"
  std::cout <<"#r val dfunc ( r_0="<<switchingFunction.description()<<")\n";
  double x=lowerLimit;
  while(x < upperLimit) {
    double dfunc;
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
