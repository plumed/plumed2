/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.0.

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

#include "FunctionVessel.h"
#include "ActionWithDistribution.h"

namespace PLMD {

class min : public SumVessel {
private:
  double beta;
public:
  static void reserveKeyword( Keywords& keys );
  min( const VesselOptions& da );
  double compute( const unsigned& i, const double& val, double& df );
  double final_computations( const unsigned& i, const double& valin, double& df );
  void printKeywords();
};

PLUMED_REGISTER_VESSEL(min,"MIN")

void min::reserveKeyword( Keywords& keys ){
  keys.reserve("optional","MIN","calculate the minimum value and store it in a value called min. "
                                "To make this quantity continuous the minimum is calculated using "
                                "\\f$ \\textrm{min} = \\frac{\\beta}{ \\log \\sum_i \\exp\\left( \\frac{\\beta}{s_i} \\right) } \\f$ "
                                "The value of \\f$\\beta\\f$ in this function is specified using (BETA=\\f$\\beta\\f$)");
}

min::min( const VesselOptions& da ) :
SumVessel(da)
{
  std::vector<std::string> data=Tools::getWords(da.parameters);
  if( data.size()!=1 ){ error("There should only be a value for beta in the input to MIN"); return; }
  bool found_beta=Tools::parse(data,"BETA",beta);
  if (!found_beta){ error("No value for BETA specified in call to MIN"); return; } 

  addOutput("min");
  log.printf("  value %s.min contains the minimum value. beta is equal to %f\n",(getAction()->getLabel()).c_str(),beta);
}

void min::printKeywords(){
  Keywords mkeys; 
  mkeys.add("compulsory","BETA","the value of beta for the equation in the manual");
  mkeys.print(log);
}

double min::compute( const unsigned& i, const double& val, double& df ){
  double f; f=exp( beta/val ); df=f/(val*val);
  return f;
}

double min::final_computations( const unsigned& i, const double& valin, double& df ){
  double dist; dist=beta/std::log( valin ); df=dist*dist/valin;
  return dist;
}

}
