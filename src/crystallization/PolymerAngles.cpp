/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2017-2019 The plumed team
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

//+PLUMEDOC MCOLVARF POLYMER_ANGLES
/*
Calculate a function to investigate the relative orientations of polymer angles

This CV takes the vectors calculated by a \ref PLANES action as input and computes the following function
of the relative angles, \f$\theta\f$, between the vectors that are normal to the pairs of input vectors:

\f[
s = \frac{ 3 \cos^2 \theta - 1 }{ 2 }
\f]

This average of this quantity over all the vectors in the first coordination sphere around each of the PLANES specified
is then calculated.

\par Examples

The example below calculates a set of vectors using the \ref PLANES action.  The average number for the function \f$s\f$
defined above is then computed over the first coordination sphere of each of the centers of mass of the molecules that were
used to define the planes.  Finally the average of these quantities is computed an printed to a file.

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
LABEL=m3
... PLANES

s3: POLYMER_ANGLES SPECIES=m3 LOWMEM SWITCH={RATIONAL R_0=0.6} MEAN
PRINT ARG=s3.mean FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace crystallization {

class PolymerAngles : public OrientationSphere {
public:
  static void registerKeywords( Keywords& keys );
  explicit PolymerAngles(const ActionOptions& ao);
  double computeVectorFunction( const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
                                Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const ;
};

PLUMED_REGISTER_ACTION(PolymerAngles,"POLYMER_ANGLES")

void PolymerAngles::registerKeywords( Keywords& keys ) {
  OrientationSphere::registerKeywords(keys);
}

PolymerAngles::PolymerAngles(const ActionOptions& ao):
  Action(ao),
  OrientationSphere(ao)
{
  if( mybasemulticolvars.size()==0 ) error("SMAC must take multicolvar as input");
  for(unsigned i=0; i<mybasemulticolvars.size(); ++i) {
    if( (mybasemulticolvars[i]->getNumberOfQuantities()-2)%3!=0 ) error("POLYMER_ANGLES is only possible with three dimensional vectors");
  }
}

double PolymerAngles::computeVectorFunction( const Vector& conn, const std::vector<double>& vec1, const std::vector<double>& vec2,
    Vector& dconn, std::vector<double>& dvec1, std::vector<double>& dvec2 ) const {

  plumed_assert( (vec1.size()-2)==3 );
  double dot = 0; for(unsigned k=0; k<3; ++k) dot += vec1[2+k]*vec2[2+k];
  double ans = 1.5*dot*dot - 0.5; for(unsigned k=0; k<3; ++k) { dvec1[2+k]=3*dot*vec2[2+k]; dvec2[2+k]=3*dot*vec1[2+k]; }
  return ans;
}

}
}
