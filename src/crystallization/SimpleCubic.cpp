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
#include "CubicHarmonicBase.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace crystallization {

//+PLUMEDOC MCOLVAR SIMPLECUBIC
/*
Calculate whether or not the coordination spheres of atoms are arranged as they would be in a simple cubic structure.

We can measure how similar the environment around atom \f$i\f$ is to a simple cubic structure is by evaluating
the following quantity:

\f[
s_i = \frac{ \sum_{i \ne j} \sigma(r_{ij}) \left[ \frac{ x_{ij}^4 + y_{ij}^4 + z_{ij}^4 }{r_{ij}^4} \right] }{ \sum_{i \ne j} \sigma(r_{ij}) }
\f]

In this expression \f$x_{ij}\f$, \f$y_{ij}\f$ and \f$z_{ij}\f$ are the \f$x\f$, \f$y\f$ and \f$z\f$ components of the vector connecting atom \f$i\f$ to
atom \f$j\f$ and \f$r_{ij}\f$ is the magnitude of this vector.  \f$\sigma(r_{ij})\f$ is a \ref switchingfunction that acts on the distance between atom \f$i\f$ and atom \f$j\f$ and its inclusion in the numerator and the denominator of the above expression as well as the fact that we are summing
over all of the other atoms in the system ensures that we are calculating an average
of the function of \f$x_{ij}\f$, \f$y_{ij}\f$ and \f$z_{ij}\f$ for the atoms in the first coordination sphere around atom \f$i\f$.
This quantity is once again a multicolvar so you can compute it for multiple atoms using a single PLUMED action and then compute
the average value for the atoms in your system, the number of atoms that have an \f$s_i\f$ value that is more that some target and
so on.  Notice also that you can rotate the reference frame if you are using a non-standard unit cell.


\par Examples

The following input tells plumed to calculate the simple cubic parameter for the atoms 1-100 with themselves.
The mean value is then calculated.
\plumedfile
SIMPLECUBIC SPECIES=1-100 R_0=1.0 MEAN
\endplumedfile

The following input tells plumed to look at the ways atoms 1-100 are within 3.0 are arranged about atoms
from 101-110.  The number of simple cubic parameters that are greater than 0.8 is then output
\plumedfile
SIMPLECUBIC SPECIESA=101-110 SPECIESB=1-100 R_0=3.0 MORE_THAN={RATIONAL R_0=0.8 NN=6 MM=12 D_0=0}
\endplumedfile

*/
//+ENDPLUMEDOC


class SimpleCubic : public CubicHarmonicBase {
public:
  static void registerKeywords( Keywords& keys );
  explicit SimpleCubic(const ActionOptions&);
  double calculateCubicHarmonic( const Vector& distance, const double& d2, Vector& myder ) const ;
};

PLUMED_REGISTER_ACTION(SimpleCubic,"SIMPLECUBIC")

void SimpleCubic::registerKeywords( Keywords& keys ) {
  CubicHarmonicBase::registerKeywords( keys );
}

SimpleCubic::SimpleCubic(const ActionOptions&ao):
  Action(ao),
  CubicHarmonicBase(ao)
{
  checkRead();
}

double SimpleCubic::calculateCubicHarmonic( const Vector& distance, const double& d2, Vector& myder ) const {
  double x2 = distance[0]*distance[0];
  double x3 = distance[0]*x2;
  double x4 = distance[0]*x3;

  double y2 = distance[1]*distance[1];
  double y3 = distance[1]*y2;
  double y4 = distance[1]*y3;

  double z2 = distance[2]*distance[2];
  double z3 = distance[2]*z2;
  double z4 = distance[2]*z3;

  double r4 = pow( d2, 2 );
  double tmp = ( x4 + y4 + z4 ) / r4;

  double t1=(x2+y2+z2), t2=t1*t1, t3=(x4+y4+z4)/(t1*t2);
  myder[0] = 4*x3/t2-4*distance[0]*t3;
  myder[1] = 4*y3/t2-4*distance[1]*t3;
  myder[2] = 4*z3/t2-4*distance[2]*t3;
  return tmp;
}

}
}

