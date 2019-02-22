/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2018 The VES code team
   (see the PEOPLE-VES file at the root of this folder for a list of names)

   See http://www.ves-code.org for more information.

   This file is part of VES code module.

   The VES code module is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   The VES code module is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with the VES code module.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "FermiSwitchingFunction.h"

#include "tools/Tools.h"
#include "tools/Keywords.h"

#include <vector>
#include <limits>


namespace PLMD {
namespace ves {


void FermiSwitchingFunction::registerKeywords( Keywords& keys ) {
  keys.add("compulsory","R_0","the value of R_0 in the switching function");
  keys.add("compulsory","FERMI_LAMBDA","the value of lambda in the Fermi-type switching function (only needed for TYPE=FERMI).");
  keys.add("optional","FERMI_EXP_MAX","only needed for TYPE=FERMI");
}

void FermiSwitchingFunction::set(const std::string& definition,std::string& errormsg) {
  std::vector<std::string> data=Tools::getWords(definition);
  if( data.size()<1 ) {
    errormsg="missing all input for switching function";
    return;
  }
  std::string name=data[0];
  data.erase(data.begin());
  if(name!="FERMI") {errormsg="only FERMI is supported";}
  type=fermi;
  //
  bool found_r0=Tools::parse(data,"R_0",r0_);
  if(!found_r0) {errormsg="R_0 is required";}

  //
  fermi_exp_max_=std::numeric_limits<double>::max();
  Tools::parse(data,"FERMI_EXP_MAX",fermi_exp_max_);
  //
  bool found_lambda=Tools::parse(data,"FERMI_LAMBDA",fermi_lambda_);
  if(!found_lambda) {errormsg="FERMI_LAMBDA is required for FERMI";}
  if( !data.empty() ) {
    errormsg="found the following rogue keywords in switching function input : ";
    for(unsigned i=0; i<data.size(); ++i) errormsg = errormsg + data[i] + " ";
  }
  init=true;
  if(errormsg.size()>0) {init=false;}
}

std::string FermiSwitchingFunction::description() const {
  std::ostringstream ostr;
  ostr<<1./invr0_<<".  Using ";
  if(type==fermi) {
    ostr<< "fermi switching function with parameter";
    ostr<< " lambda="<<fermi_lambda_;
  }
  else {
    plumed_merror("Unknown switching function type");
  }
  return ostr.str();
}


double FermiSwitchingFunction::calculate(double distance, double& dfunc) const {
  plumed_massert(init,"you are trying to use an unset FermiSwitchingFunction");
  double rdist=fermi_lambda_*(distance-r0_);
  if(rdist >= fermi_exp_max_) {rdist = fermi_exp_max_;}
  double result = 1.0/(1.0+exp(rdist));
  dfunc=-fermi_lambda_*exp(rdist)*result*result;
  // this is because calculate() sets dfunc to the derivative divided times the distance.
  // (I think this is misleading and I would like to modify it - GB)
  // dfunc/=distance;
  //
  return result;
}


FermiSwitchingFunction::FermiSwitchingFunction():
  init(false),
  type(fermi),
  r0_(0.0),
  invr0_(0.0),
  fermi_lambda_(1.0),
  fermi_exp_max_(100.0)
{
}

FermiSwitchingFunction::FermiSwitchingFunction(const FermiSwitchingFunction&sf):
  init(sf.init),
  type(sf.type),
  r0_(sf.r0_),
  invr0_(sf.invr0_),
  fermi_lambda_(sf.fermi_lambda_),
  fermi_exp_max_(sf.fermi_exp_max_)
{
}

void FermiSwitchingFunction::set(const double r0, const double fermi_lambda, const double fermi_exp_max) {
  init=true;
  type=fermi;
  r0_=r0;
  fermi_lambda_=fermi_lambda;
  if(fermi_exp_max>0.0) {
    fermi_exp_max_=fermi_exp_max;
  }
  else {
    fermi_exp_max_=100.0;
  }

}

double FermiSwitchingFunction::get_r0() const {
  return r0_;
}


FermiSwitchingFunction::~FermiSwitchingFunction() {
}


}
}
