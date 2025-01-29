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
#include "FunctionTemplateBase.h"
#include "FunctionShortcut.h"
#include "FunctionOfScalar.h"

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION PIECEWISE
/*
Compute a piece wise straight line through its arguments that passes through a set of ordered control points.

For variables less than the first
(greater than the last) point, the value of the first (last) point is used.

\f[
\frac{y_{i+1}-y_i}{x_{i+1}-x_i}(s-x_i)+y_i ;  if x_i<s<x_{i+1}
\f]
\f[
y_N ; if x>x_{N-1}
\f]
\f[
y_1 ; if x<x_0
\f]

Control points are passed using the POINT0=... POINT1=... syntax as in the example below

If one argument is supplied, it results in a scalar quantity.
If multiple arguments are supplied, it results
in a vector of values. Each value will be named as the name of the original
argument with suffix _pfunc.

\par Examples

\plumedfile
dist1: DISTANCE ATOMS=1,10
dist2: DISTANCE ATOMS=2,11

pw: PIECEWISE POINT0=1,10 POINT1=2,PI POINT2=3,10 ARG=dist1
ppww: PIECEWISE POINT0=1,10 POINT1=2,PI POINT2=3,10 ARG=dist1,dist2
PRINT ARG=pw,ppww.dist1_pfunc,ppww.dist2_pfunc
\endplumedfile


*/
//+ENDPLUMEDOC

//+PLUMEDOC FUNCTION PIECEWISE_SCALAR
/*
Compute a piece wise straight line through its arguments that passes through a set of ordered control points.

\par Examples

*/
//+ENDPLUMEDOC

class Piecewise : public FunctionTemplateBase {
  std::vector<std::pair<double,double> > points;
public:
  void registerKeywords(Keywords& keys) override;
  void read( ActionWithArguments* action ) override;
  void setPeriodicityForOutputs( ActionWithValue* action ) override;
  void calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const override;
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

void Piecewise::read( ActionWithArguments* action ) {
  for(int i=0;; i++) {
    std::vector<double> pp;
    if(!action->parseNumberedVector("POINT",i,pp) ) {
      break;
    }
    if(pp.size()!=2) {
      action->error("points should be in x,y format");
    }
    points.push_back(std::pair<double,double>(pp[0],pp[1]));
    if(i>0 && points[i].first<=points[i-1].first) {
      action->error("points abscissas should be monotonously increasing");
    }
  }

  for(unsigned i=0; i<action->getNumberOfArguments(); i++) {
    if(action->getPntrToArgument(i)->isPeriodic()) {
      action->error("Cannot use PIECEWISE on periodic arguments");
    }
  }
  action->log.printf("  on points:");
  for(unsigned i=0; i<points.size(); i++) {
    action->log.printf("   (%f,%f)",points[i].first,points[i].second);
  }
  action->log.printf("\n");
}

void Piecewise::setPeriodicityForOutputs( ActionWithValue* action ) {
  for(unsigned i=0; i<action->getNumberOfComponents(); ++i) {
    action->copyOutput(i)->setNotPeriodic();
  }
}

void Piecewise::calc( const ActionWithArguments* action, const std::vector<double>& args, std::vector<double>& vals, Matrix<double>& derivatives ) const {
  for(unsigned i=0; i<args.size(); i++) {
    unsigned p=0;
    for(; p<points.size(); p++) {
      if(args[i]<points[p].first) {
        break;
      }
    }
    if(p==0) {
      vals[i]=points[0].second;
      derivatives(i,i)=0.0;
    } else if(p==points.size()) {
      vals[i]=points[points.size()-1].second;
      derivatives(i,i)=0.0;
    } else {
      double m=(points[p].second-points[p-1].second) / (points[p].first-points[p-1].first);
      vals[i]=m*(args[i]-points[p-1].first)+points[p-1].second;
      derivatives(i,i)=m;
    }
  }
}

}
}


