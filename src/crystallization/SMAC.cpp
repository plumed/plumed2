/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2016 The plumed team
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
#include "OrientationSphere.h"
#include "core/ActionRegister.h"
#include "tools/KernelFunctions.h"

//+PLUMEDOC MCOLVARF SMAC
/*
Calculate the SMAC collective variable discussed in \cite smac-paper

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class SMAC : public OrientationSphere {
private:
  std::vector<KernelFunctions> kernels;
  SwitchingFunction coord_switch;
public:
  static void registerKeywords( Keywords& keys ); 
  explicit SMAC(const ActionOptions& ao); 
  double transformDotProduct( const double& dot, double& df ) const ; 
  double calculateCoordinationPrefactor( const double& coord, double& df ) const ;
};

PLUMED_REGISTER_ACTION(SMAC,"SMAC")

void SMAC::registerKeywords( Keywords& keys ){
  OrientationSphere::registerKeywords(keys);
  keys.add("numbered","KERNEL","The kernels used in the function of the angle");
  keys.add("compulsory","SWITCH_COORD","This keyword is used to define the coordination switching function.");
  keys.reset_style("KERNEL","compulsory"); 
}

SMAC::SMAC(const ActionOptions& ao):
Action(ao),
OrientationSphere(ao)
{
   std::string kernelinpt;
   for(int i=1;;i++){
      if( !parseNumbered("KERNEL",i,kernelinpt) ) break;
      KernelFunctions mykernel( kernelinpt, false );
      kernels.push_back( mykernel ); 
   }
   if( kernels.size()==0 ) error("no kernels defined");

   std::string sw, errors; parse("SWITCH_COORD",sw);
   if(sw.length()==0) error("SWITCH_COORD keyword is missing");
   coord_switch.set(sw,errors);
   if(errors.length()>0) error("the following errors were found in input to SWITCH_COORD : " + errors );

}

double SMAC::transformDotProduct( const double& dot, double& df ) const {
  std::vector<Value*> pos; pos.push_back( new Value() ); std::vector<double> deriv(1);
  pos[0]->setNotPeriodic(); pos[0]->set( acos( dot ) ); 
  double ans=0; df=0; double dcos=-1./sqrt( 1. - dot*dot );
  for(unsigned i=0;i<kernels.size();++i){
      ans += kernels[i].evaluate( pos, deriv );
      df += deriv[0]*dcos;
  }
  delete pos[0];
  return ans;
}

double SMAC::calculateCoordinationPrefactor( const double& coord, double& df ) const {
  double f=coord_switch.calculate( coord, df ); df*=coord; return f;
}

}
}
