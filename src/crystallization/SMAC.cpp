/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

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

namespace PLMD {
namespace crystallization {

class SMAC : public OrientationSphere {
private:
  std::vector<KernelFunctions> kernels;
  std::vector<double> deriv;
  std::vector<Value*> pos;
public:
  static void registerKeywords( Keywords& keys ); 
  SMAC(const ActionOptions& ao); 
  ~SMAC();
  double transformDotProduct( const double& dot, double& df ); 
};

PLUMED_REGISTER_ACTION(SMAC,"SMAC")

void SMAC::registerKeywords( Keywords& keys ){
  OrientationSphere::registerKeywords(keys);
  keys.add("numbered","KERNEL","The kernels used in the function of the angle");
  keys.reset_style("KERNEL","compulsory"); 
}

SMAC::SMAC(const ActionOptions& ao):
Action(ao),
OrientationSphere(ao),
deriv(1)
{
   std::string kernelinpt;
   for(int i=1;;i++){
      if( !parseNumbered("KERNEL",i,kernelinpt) ) break;
      KernelFunctions mykernel( kernelinpt, false );
      kernels.push_back( mykernel ); 
   }
   if( kernels.size()==0 ) error("no kernels defined");

   pos.push_back( new Value() ); 
   pos[0]->setNotPeriodic();
}

SMAC::~SMAC(){
   delete pos[0];
}

double SMAC::transformDotProduct( const double& dot, double& df ){
  double ans=0; df=0; pos[0]->set( acos( dot ) ); double dcos=-1./sqrt( 1. - dot*dot );
  for(unsigned i=0;i<kernels.size();++i){
      ans += kernels[i].evaluate( pos, deriv );
      df += deriv[0]*dcos;
  }
  return ans;
}

}
}
