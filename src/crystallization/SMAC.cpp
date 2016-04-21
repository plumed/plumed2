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
#include "tools/Torsion.h"
#include "tools/KernelFunctions.h"

//+PLUMEDOC MCOLVARF SMAC
/*
Calculate the SMAC collective variable discussed in \cite smac-paper


\bug Contribution to virial is wrong

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
  double computeVectorFunction( const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2, 
                                Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const ;
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
   if( mybasemulticolvars.size()==0 ) error("SMAC must take multicolvar as input");
   for(unsigned i=0;i<mybasemulticolvars.size();++i){
       if( (mybasemulticolvars[i]->getNumberOfQuantities()-2)%3!=0 ) error("SMAC is only possible with three dimensional vectors");
   }

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

double SMAC::computeVectorFunction( const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2, 
                                    Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const {

  unsigned nvectors = ( vec1.size() - 2 ) / 3; plumed_assert( (vec1.size()-2)%3==0 );
  std::vector<Vector> dv1(nvectors), dv2(nvectors), tdconn(nvectors); Torsion t; std::vector<Vector> v1(nvectors), v2(nvectors);
  std::vector<Value*> pos; for(unsigned i=0;i<nvectors;++i){ pos.push_back( new Value() ); pos[i]->setDomain( "-pi", "pi" ); }

  for(unsigned j=0;j<nvectors;++j){
      for(unsigned k=0;k<3;++k){
          v1[j][k]=vec1[2+3*j+k]; v2[j][k]=vec2[2+3*j+k];
      }
      double angle = t.compute( v1[j], conn, v2[j], dv1[j], tdconn[j], dv2[j] );
      pos[j]->set( angle );
  }

  double ans=0; std::vector<double> deriv( nvectors ), df( nvectors, 0 ); 
  for(unsigned i=0;i<kernels.size();++i){
      ans += kernels[i].evaluate( pos, deriv );
      for(unsigned j=0;j<nvectors;++j) df[j] += deriv[j];
  }
  dconn.zero(); for(unsigned j=0;j<nvectors;++j) dconn += df[j]*tdconn[j];
  for(unsigned j=0;j<nvectors;++j){
      for(unsigned k=0;k<3;++k){ dvec1[2+3*j+k]=df[j]*dv1[j][k]; dvec2[2+3*j+k]=df[j]*dv2[j][k]; }
      delete pos[j];
  }
  return ans;
}

double SMAC::calculateCoordinationPrefactor( const double& coord, double& df ) const {
  double f=1-coord_switch.calculate( coord, df ); df*=-coord; return f;
}

}
}
