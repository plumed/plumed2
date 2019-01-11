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

class Max : public FunctionVessel {
private:
  double beta;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit Max( const VesselOptions& da );
  std::string value_descriptor();
  double calcTransform( const double& val, double& dv ) const ;
  double finalTransform( const double& val, double& dv );
};

PLUMED_REGISTER_VESSEL(Max,"MAX")

void Max::registerKeywords( Keywords& keys ) {
  FunctionVessel::registerKeywords( keys );
  keys.add("compulsory","BETA","the value of beta for the equation in the manual");
}

void Max::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","MAX","calculate the maximum value. "
               "To make this quantity continuous the maximum is calculated using "
               "\\f$ \\textrm{max} = \\beta \\log \\sum_i \\exp\\left( \\frac{s_i}{\\beta}\\right) \\f$ "
               "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
  keys.addOutputComponent("max","MAX","the maximum value. This is calculated using the formula described in the description of the "
                          "keyword so as to make it continuous.");

}

Max::Max( const VesselOptions& da ) :
  FunctionVessel(da)
{
  if( getAction()->isPeriodic() ) error("max is not a meaningful option for periodic variables");
  parse("BETA",beta);

  if( diffweight ) error("can't calculate max if weight is differentiable");
}

std::string Max::value_descriptor() {
  std::string str_beta; Tools::convert( beta, str_beta );
  return "the maximum value. Beta is equal to " + str_beta;
}

double Max::calcTransform( const double& val, double& dv ) const {
  double f = exp(val/beta); dv=f/beta; return f;
}

double Max::finalTransform( const double& val, double& dv ) {
  double dist=beta*std::log( val );
  dv = beta/val; return dist;
}

}
}
