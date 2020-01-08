/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2020 The plumed team
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

using namespace std;

namespace PLMD {
namespace function {

//+PLUMEDOC FUNCTION STATS
/*
Calculates statistical properties of a set of collective variables with respect to a set of reference values.

In particular it calculates and stores as components the sum of the squared deviations, the correlation, the
slope and the intercept of a linear fit.

The reference values can be either provided as values using PARAMETERS or using value without derivatives
from other actions using PARARG (for example using experimental values from collective variables such as
\ref CS2BACKBONE, \ref RDC, \ref NOE, \ref PRE).

\par Examples

The following input tells plumed to print the distance between three couple of atoms
and compare them with three reference distances.

\plumedfile
d1: DISTANCE ATOMS=10,50
d2: DISTANCE ATOMS=1,100
d3: DISTANCE ATOMS=45,75
st: STATS ARG=d1,d2,d3 PARAMETERS=1.5,4.0,2.0
PRINT ARG=d1,d2,d3,st.*
\endplumedfile

*/
//+ENDPLUMEDOC


class Stats :
  public Function
{
  std::vector<double> parameters;
  bool sqdonly;
  bool components;
  bool upperd;
public:
  explicit Stats(const ActionOptions&);
  void calculate() override;
  static void registerKeywords(Keywords& keys);
};


PLUMED_REGISTER_ACTION(Stats,"STATS")

void Stats::registerKeywords(Keywords& keys) {
  Function::registerKeywords(keys);
  keys.use("ARG");
  keys.add("optional","PARARG","the input for this action is the scalar output from one or more other actions without derivatives.");
  keys.add("optional","PARAMETERS","the parameters of the arguments in your function");
  keys.addFlag("SQDEVSUM",false,"calculates only SQDEVSUM");
  keys.addFlag("SQDEV",false,"calculates and store the SQDEV as components");
  keys.addFlag("UPPERDISTS",false,"calculates and store the SQDEV as components");
  keys.addOutputComponent("sqdevsum","default","the sum of the squared deviations between arguments and parameters");
  keys.addOutputComponent("corr","default","the correlation between arguments and parameters");
  keys.addOutputComponent("slope","default","the slope of a linear fit between arguments and parameters");
  keys.addOutputComponent("intercept","default","the intercept of a linear fit between arguments and parameters");
  keys.addOutputComponent("sqd","SQDEV","the squared deviations between arguments and parameters");
}

Stats::Stats(const ActionOptions&ao):
  Action(ao),
  Function(ao),
  sqdonly(false),
  components(false),
  upperd(false)
{
  parseVector("PARAMETERS",parameters);
  if(parameters.size()!=static_cast<unsigned>(getNumberOfArguments())&&!parameters.empty())
    error("Size of PARAMETERS array should be either 0 or the same as of the number of arguments in ARG1");

  vector<Value*> arg2;
  parseArgumentList("PARARG",arg2);

  if(!arg2.empty()) {
    if(parameters.size()>0) error("It is not possible to use PARARG and PARAMETERS together");
    if(arg2.size()!=getNumberOfArguments()) error("Size of PARARG array should be the same as number for arguments in ARG");
    for(unsigned i=0; i<arg2.size(); i++) {
      parameters.push_back(arg2[i]->get());
      if(arg2[i]->hasDerivatives()==true) error("PARARG can only accept arguments without derivatives");
    }
  }

  if(parameters.size()!=getNumberOfArguments())
    error("PARARG or PARAMETERS arrays should include the same number of elements as the arguments in ARG");

  if(getNumberOfArguments()<2) error("STATS need at least two arguments to be used");

  parseFlag("SQDEVSUM",sqdonly);
  parseFlag("SQDEV",components);
  parseFlag("UPPERDISTS",upperd);

  if(sqdonly&&components) error("You cannot used SQDEVSUM and SQDEV at the sametime");

  if(components) sqdonly = true;

  if(!arg2.empty()) log.printf("  using %zu parameters from inactive actions:", arg2.size());
  else              log.printf("  using %zu parameters:", arg2.size());
  for(unsigned i=0; i<parameters.size(); i++) log.printf(" %f",parameters[i]);
  log.printf("\n");

  if(sqdonly) {
    if(components) {
      for(unsigned i=0; i<parameters.size(); i++) {
        std::string num; Tools::convert(i,num);
        addComponentWithDerivatives("sqd-"+num);
        componentIsNotPeriodic("sqd-"+num);
      }
    } else {
      addComponentWithDerivatives("sqdevsum");
      componentIsNotPeriodic("sqdevsum");
    }
  } else {
    addComponentWithDerivatives("sqdevsum");
    componentIsNotPeriodic("sqdevsum");
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
    Value* val;
    if(!components) val=getPntrToComponent("sqdevsum");
    for(unsigned i=0; i<parameters.size(); ++i) {
      double dev = getArgument(i)-parameters[i];
      if(upperd&&dev<0) dev=0.;
      if(components) {
        val=getPntrToComponent(i);
        val->set(dev*dev);
      } else {
        nsqd += dev*dev;
      }
      setDerivative(val,i,2.*dev);
    }
    if(!components) val->set(nsqd);

  } else {

    double scx=0., scx2=0., scy=0., scy2=0., scxy=0.;

    for(unsigned i=0; i<parameters.size(); ++i) {
      const double tmpx=getArgument(i);
      const double tmpy=parameters[i];
      scx  += tmpx;
      scx2 += tmpx*tmpx;
      scy  += tmpy;
      scy2 += tmpy*tmpy;
      scxy += tmpx*tmpy;
    }

    const double ns = parameters.size();

    const double num = ns*scxy - scx*scy;
    const double idev2x = 1./(ns*scx2-scx*scx);
    const double idevx = sqrt(idev2x);
    const double idevy = 1./sqrt(ns*scy2-scy*scy);

    /* sd */
    const double nsqd = scx2 + scy2 - 2.*scxy;
    /* correlation */
    const double correlation = num * idevx * idevy;
    /* slope and intercept */
    const double slope = num * idev2x;
    const double inter = (scy - slope * scx)/ns;

    Value* valuea=getPntrToComponent("sqdevsum");
    Value* valueb=getPntrToComponent("corr");
    Value* valuec=getPntrToComponent("slope");
    Value* valued=getPntrToComponent("intercept");

    valuea->set(nsqd);
    valueb->set(correlation);
    valuec->set(slope);
    valued->set(inter);

    /* derivatives */
    for(unsigned i=0; i<parameters.size(); ++i) {
      const double common_d1 = (ns*parameters[i]-scy)*idevx;
      const double common_d2 = num*(ns*getArgument(i)-scx)*idev2x*idevx;
      const double common_d3 = common_d1 - common_d2;

      /* sqdevsum */
      const double sq_der = 2.*(getArgument(i)-parameters[i]);
      /* correlation */
      const double co_der = common_d3*idevy;
      /* slope */
      const double sl_der = (common_d1-2.*common_d2)*idevx;
      /* intercept */
      const double int_der = -(slope+ scx*sl_der)/ns;

      setDerivative(valuea,i,sq_der);
      setDerivative(valueb,i,co_der);
      setDerivative(valuec,i,sl_der);
      setDerivative(valued,i,int_der);
    }

  }
}

}
}


