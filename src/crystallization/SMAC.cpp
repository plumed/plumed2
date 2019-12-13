/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2019 The plumed team
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
Calculate a variant on the SMAC collective variable discussed in \cite smac-paper

The SMAC collective variable can be used to study the formation of molecular solids
from either the melt or from solution.  The idea behind this variable is that what
differentiates a molecular solid from a molecular liquid is an alignment of
internal vectors in neighboring molecules.  In other words, the relative orientation
of neighboring molecules is no longer random as it is in a liquid.  In a solid particular
torsional angles between molecules are preferred.  As such this CV calculates the following
average:

\f[
s_i = \frac{ \left\{ 1 - \psi\left[ \sum_{j \ne i} \sigma(r_{ij}) \right] \right\} \sum_{j \ne i} \sigma(r_{ij}) \sum_n K_n(\theta_{ij}) }{ \sum_{j \ne i} \sigma(r_{ij}) }
\f]

In this expression \f$r_{ij}\f$ is the distance between molecule \f$i\f$ and molecule \f$j\f$ and \f$\sigma(r_{ij})\f$ is a
\ref switchingfunction that acts on this distance.  By including this switching function in the second summation in the
numerator and in the denominator we are thus ensuring that we calculate an average over the molecules in the first coordination
sphere of molecule \f$i\f$.  All molecules in higher coordination sphere will essentially contribute zero to the sums in the
above expression because their \f$\sigma(r_{ij})\f$ will be very small.  \f$\psi\f$ is also a switching function.  The term
including \f$\psi\f$ in the numerator is there to ensure that only those molecules that are attached to a reasonably large
number of molecules.  It is important to include this "more than" switching function when you are simulating nucleation
from solution with this CV.  Lastly, the $K_n functions are \ref kernelfunctions that take the torsion angle, \f$\theta_{ij}\f$, between the
internal orientation vectors for molecules \f$i\f$ and \f$j\f$ as input.  These kernel functions should be set so that they are
equal to one when the relative orientation of the molecules are as they are in the solid and equal to zero otherwise.
The final \f$s_i\f$ quantity thus measures whether (on average) the molecules in the first coordination sphere around molecule \f$i\f$
are oriented as they would be in the solid.  Furthermore, this Action is a multicolvar so you can calculate the \f$s_i\f$ values
for all the molecules in your system simultaneously and then determine the average, the number less than and so on.

\par Examples

In the example below the orientation of the molecules in the system is determined by calculating the
vector that connects a pair of atoms.  SMAC is then used to determine whether the molecules are sitting
in a solid or liquid like environment.  We can determine whether the environment is solid or liquid like because in the solid the torsional angle between
the bond vectors on adjacent molecules is close to 0 or \f$\pi\f$.  The final quantity that is output to the colvar
file measures the number of molecules that have a SMAC parameter that is greater than 0.7.  N.B. By using
the indices of three atoms for each of the MOL keywords below we are telling PLUMED to use the first two
numbers to determine the orientation of the molecule that will ultimately be used when calculating the \f$\theta_{ij}\f$
terms in the formula above.  The atom with the third index meanwhile is used when we calculate \f$r_{ij}\f$.

\plumedfile
MOLECULES ...
MOL1=9,10,9
MOL2=89,90,89
MOL3=473,474,473
MOL4=1161,1162,1161
MOL5=1521,1522,1521
MOL6=1593,1594,1593
MOL7=1601,1602,1601
MOL8=2201,2202,2201
LABEL=m3
... MOLECULES

SMAC ...
   SPECIES=m3 LOWMEM
   KERNEL1={GAUSSIAN CENTER=0 SIGMA=0.480} KERNEL2={GAUSSIAN CENTER=pi SIGMA=0.480}
   SWITCH={RATIONAL R_0=0.6} MORE_THAN={RATIONAL R_0=0.7} SWITCH_COORD={EXP R_0=4}
   LABEL=s2
... SMAC

PRINT ARG=s2.* FILE=colvar
\endplumedfile

This second example works in a way that is very similar to the previous command.  Now, however,
the orientation of the molecules is determined by finding the plane that contains the positions
of three atoms.

\plumedfile
PLANES ...
MOL1=9,10,11
MOL2=89,90,91
MOL3=473,474,475
MOL4=1161,1162,1163
MOL5=1521,1522,1523
MOL6=1593,1594,1595
MOL7=1601,1602,1603
MOL8=2201,2202,2203
VMEAN
LABEL=m3
... PLANES

SMAC ...
   SPECIES=m3 LOWMEM
   KERNEL1={GAUSSIAN CENTER=0 SIGMA=0.480} KERNEL2={GAUSSIAN CENTER=pi SIGMA=0.480}
   SWITCH={RATIONAL R_0=0.6} MORE_THAN={RATIONAL R_0=0.7} SWITCH_COORD={EXP R_0=3.0}
   LABEL=s2
... SMAC

PRINT ARG=s2.* FILE=colvar
\endplumedfile

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
                                Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const override;
  double calculateCoordinationPrefactor( const double& coord, double& df ) const override;
};

PLUMED_REGISTER_ACTION(SMAC,"SMAC")

void SMAC::registerKeywords( Keywords& keys ) {
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
  for(unsigned i=0; i<mybasemulticolvars.size(); ++i) {
    if( (mybasemulticolvars[i]->getNumberOfQuantities()-2)%3!=0 ) error("SMAC is only possible with three dimensional vectors");
  }

  std::string kernelinpt;
  for(int i=1;; i++) {
    if( !parseNumbered("KERNEL",i,kernelinpt) ) break;
    KernelFunctions mykernel( kernelinpt );
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
  std::vector<std::unique_ptr<Value>> pos;
  for(unsigned i=0; i<nvectors; ++i) { pos.emplace_back( new Value() ); pos[i]->setDomain( "-pi", "pi" ); }

  for(unsigned j=0; j<nvectors; ++j) {
    for(unsigned k=0; k<3; ++k) {
      v1[j][k]=vec1[2+3*j+k]; v2[j][k]=vec2[2+3*j+k];
    }
    double angle = t.compute( v1[j], conn, v2[j], dv1[j], tdconn[j], dv2[j] );
    pos[j]->set( angle );
  }

  auto pos_ptr=Tools::unique2raw(pos);

  double ans=0; std::vector<double> deriv( nvectors ), df( nvectors, 0 );
  for(unsigned i=0; i<kernels.size(); ++i) {
    ans += kernels[i].evaluate( pos_ptr, deriv );
    for(unsigned j=0; j<nvectors; ++j) df[j] += deriv[j];
  }
  dconn.zero(); for(unsigned j=0; j<nvectors; ++j) dconn += df[j]*tdconn[j];
  for(unsigned j=0; j<nvectors; ++j) {
    for(unsigned k=0; k<3; ++k) { dvec1[2+3*j+k]=df[j]*dv1[j][k]; dvec2[2+3*j+k]=df[j]*dv2[j][k]; }
  }
  return ans;
}

double SMAC::calculateCoordinationPrefactor( const double& coord, double& df ) const {
  double f=1-coord_switch.calculate( coord, df ); df*=-coord; return f;
}

}
}
