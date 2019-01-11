/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2015-2019 The plumed team
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

namespace PLMD {
namespace vesselbase {

class AltMin : public vesselbase::FunctionVessel {
private:
  double beta;
public:
  static void registerKeywords( Keywords& keys );
  static void reserveKeyword( Keywords& keys );
  explicit AltMin( const vesselbase::VesselOptions& da );
  std::string value_descriptor();
  double calcTransform( const double& val, double& dv ) const ;
  double finalTransform( const double& val, double& dv );
};

PLUMED_REGISTER_VESSEL(AltMin,"ALT_MIN")

void AltMin::registerKeywords( Keywords& keys ) {
  FunctionVessel::registerKeywords(keys);
  keys.add("compulsory","BETA","the value of beta for the equation in the manual");
}

void AltMin::reserveKeyword( Keywords& keys ) {
  keys.reserve("vessel","ALT_MIN","calculate the minimum value. "
               "To make this quantity continuous the minimum is calculated using "
               "\\f$ \\textrm{min} = -\\frac{1}{\\beta} \\log \\sum_i \\exp\\left( -\\beta s_i \\right)  \\f$ "
               "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$).");
  keys.addOutputComponent("altmin","ALT_MIN","the minimum value. This is calculated using the formula described in the description of the "
                          "keyword so as to make it continuous.");
}

AltMin::AltMin( const vesselbase::VesselOptions& da ):
  FunctionVessel(da)
{
  if( getAction()->isPeriodic() ) error("MIN is not a meaningful option for periodic variables");
  parse("BETA",beta); usetol=true;
}

std::string AltMin::value_descriptor() {
  std::string str_beta; Tools::convert( beta, str_beta );
  return "the minimum value. Beta is equal to " + str_beta;
}

double AltMin::calcTransform( const double& val, double& dv ) const {
  double f = exp( -beta*val ); dv = -beta*f; return f;
}

double AltMin::finalTransform( const double& val, double& dv ) {
  dv = - 1.0 /(beta*val); return -std::log( val ) / beta;
}

}
}
