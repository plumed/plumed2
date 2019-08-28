/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2019 The plumed team
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
#include "VesselRegister.h"
#include "FunctionVessel.h"
#include "ActionWithVessel.h"

namespace PLMD {
namespace vesselbase {

class Min : public FunctionVessel {
private:
  double beta;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit Min( const VesselOptions& da );
  std::string value_descriptor() override;
  double calcTransform( const double& val, double& dv ) const override;
  double finalTransform( const double& val, double& dv ) override;
};

PLUMED_REGISTER_VESSEL(Min,"MIN")

void Min::registerKeywords( Keywords& keys ) {
  FunctionVessel::registerKeywords( keys );
  keys.add("compulsory","BETA","the value of beta for the equation in the manual");
}

void Min::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","MIN","calculate the minimum value. "
               "To make this quantity continuous the minimum is calculated using "
               "\\f$ \\textrm{min} = \\frac{\\beta}{ \\log \\sum_i \\exp\\left( \\frac{\\beta}{s_i} \\right) } \\f$ "
               "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.addOutputComponent("min","MIN","the minimum value. This is calculated using the formula described in the description of the "
                          "keyword so as to make it continuous.");
}

Min::Min( const VesselOptions& da ) :
  FunctionVessel(da)
{
  if( getAction()->isPeriodic() ) error("min is not a meaningful option for periodic variables");
  parse("BETA",beta);

  if( diffweight ) error("can't calculate min if weight is differentiable");
}

std::string Min::value_descriptor() {
  std::string str_beta; Tools::convert( beta, str_beta );
  return "the minimum value. Beta is equal to " + str_beta;
}

double Min::calcTransform( const double& val, double& dv ) const {
  double f = exp(beta/val); dv=f/(val*val);
  return f;
}

double Min::finalTransform( const double& val, double& dv ) {
  double dist=beta/std::log( val );
  dv = dist*dist/val; return dist;
}

}
}
