/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015,2016 The plumed team
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
#include "Function.h"
#include "ActionRegister.h"

#include <cmath>

using namespace std;

namespace PLMD{
namespace function{

//+PLUMEDOC FUNCTION STATS
/*

\endverbatim

*/
//+ENDPLUMEDOC


class Stats :
  public Function
{
  std::vector<double> parameters;
  bool sqdonly;
public:
  explicit Stats(const ActionOptions&);
  void calculate();
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Stats,"STATS")

void Stats::registerKeywords(Keywords& keys){
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","PARARG","the input for this action is the scalar output from one or more other actions. The particular scalars that you will use "
                                  "are referenced using the label of the action. If the label appears on its own then it is assumed that the Action calculates "
                                  "a single scalar value.  The value of this scalar is thus used as the input to this new action.  If * or *.* appears the "
                                  "scalars calculated by all the proceding actions in the input file are taken.  Some actions have multi-component outputs and "
                                  "each component of the output has a specific label.  For example a \\ref DISTANCE action labelled dist may have three componets "
                                  "x, y and z.  To take just the x component you should use dist.x, if you wish to take all three components then use dist.*."
                                  "More information on the referencing of Actions can be found in the section of the manual on the PLUMED \\ref Syntax.  "
                                  "Scalar values can also be "
                                  "referenced using POSIX regular expressions as detailed in the section on \\ref Regex. To use this feature you you must compile "
                                  "PLUMED with the appropriate flag."); 
  keys.add("optional","PARAMETERS","the parameters of the arguments in your function");
  keys.addFlag("SQDEVSUM",false,"calculates only SQDEVSUM");
  keys.addOutputComponent("sqdevsum","default","the sum of the squared deviations between arguments and parameters"); 
  keys.addOutputComponent("corr","default","the correlation between arguments and parameters"); 
  keys.addOutputComponent("slope","default","the slope of a linear fit between arguments and parameters"); 
  keys.addOutputComponent("intercept","default","the intercept of a linear fit between arguments and parameters"); 
}

Stats::Stats(const ActionOptions&ao):
Action(ao),
Function(ao),
sqdonly(false)
{
  parseVector("PARAMETERS",parameters);
  if(parameters.size()!=static_cast<unsigned>(getNumberOfArguments())&&!parameters.empty())
    error("Size of PARAMETERS array should be either 0 or the same as of the number of arguments in ARG1");

  vector<Value*> arg2;
  parseArgumentList("PARARG",arg2);

  if(!arg2.empty()) {
    if(parameters.size()>0) error("It is not possible to use PARARG and PARAMETERS together");
    if(arg2.size()!=getNumberOfArguments()) error("Size of PARARG array should be the same as number for arguments in ARG");
    for(unsigned i=0;i<arg2.size();i++){
      parameters.push_back(arg2[i]->get()); 
      if(arg2[i]->hasDerivatives()==true) error("PARARG can only accept arguments without derivatives");
    }
  }

  parseFlag("SQDEVSUM",sqdonly);

  if(!arg2.empty()) log.printf("  using %zu parameters from inactive actions:", arg2.size());
  else              log.printf("  using %zu parameters:", arg2.size());
  for(unsigned i=0;i<parameters.size();i++) log.printf(" %f",parameters[i]);
  log.printf("\n");

  addComponentWithDerivatives("sqdevsum");
  componentIsNotPeriodic("sqdevsum");
  if(!sqdonly) {
    addComponentWithDerivatives("corr");
    componentIsNotPeriodic("corr");
    addComponentWithDerivatives("slope");
    componentIsNotPeriodic("slope");
    addComponentWithDerivatives("intercept");
    componentIsNotPeriodic("intercept");
  }

  checkRead();
}

void Stats::calculate()
{

  if(sqdonly) {

    double nsqd = 0.;
    Value* valuea=getPntrToComponent("sqdevsum");
    for(unsigned i=0;i<parameters.size();++i){
      double dev = getArgument(i)-parameters[i];
      setDerivative(valuea,i,2.*dev);
      nsqd += dev*dev;
    }
    valuea->set(nsqd);

  } else {

    double scx=0., scx2=0., scy=0., scy2=0., scxy=0.;
 
    for(unsigned i=0;i<parameters.size();++i){
      scx  += getArgument(i);
      scx2 += getArgument(i)*getArgument(i);
      scy  += parameters[i];
      scy2 += parameters[i]*parameters[i];
      scxy += getArgument(i)*parameters[i];
    }
  
    double ns = parameters.size();

    double num = ns*scxy - scx*scy;
    double idev2x = 1./(ns*scx2-scx*scx);
    double idevx = sqrt(idev2x);
    double idevy = 1./sqrt(ns*scy2-scy*scy);

    /* sd */
    double nsqd = scx2 + scy2 - 2.*scxy;
    /* correlation */
    double correlation = num * idevx * idevy;
    /* slope and intercept */
    double slope = num * idev2x;
    double inter = (scy - slope * scx)/ns;

    Value* valuea=getPntrToComponent("sqdevsum");
    Value* valueb=getPntrToComponent("corr");
    Value* valuec=getPntrToComponent("slope");
    Value* valued=getPntrToComponent("intercept");

    valuea->set(nsqd);
    valueb->set(correlation);
    valuec->set(slope);
    valued->set(inter);

    /* derivatives */
    for(unsigned i=0;i<parameters.size();++i){
      double common_d1 = (ns*parameters[i]-scy)*idevx;
      double common_d2 = num*(ns*getArgument(i)-scx)*idev2x*idevx;
      double common_d3 = common_d1 - common_d2;

      /* sqdevsum */
      double sq_der = 2.*(getArgument(i)-parameters[i]);
      /* correlation */
      double co_der = common_d3*idevy;
      /* slope */
      double sl_der = (common_d1-2.*common_d2)*idevx;
      /* intercept */
      double int_der = -(slope+ scx*sl_der)/ns;

      setDerivative(valuea,i,sq_der);
      setDerivative(valueb,i,co_der);
      setDerivative(valuec,i,sl_der);
      setDerivative(valued,i,int_der);
    }

  }
}

}
}


