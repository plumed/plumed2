/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014-2020 The plumed team
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
#include "tools/SwitchingFunction.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace crystallization {

//+PLUMEDOC MCOLVAR FCCUBIC
/*
Measure how similar the environment around atoms is to that found in a FCC structure.

This CV was introduced in this article \cite fcc-michele-1 and again in this article \cite fcc-michele-2
This CV essentially determines whether the environment around any given atom is similar to that found in
the FCC structure or not.  The function that is used to make this determination is as follows:

\f[
s_i = \frac{ \sum_{i \ne j} \sigma(r_{ij}) \left\{ a\left[ \frac{(x_{ij}y_{ij})^4 + (x_{ij}z_{ij})^4 + (y_{ij}z_{ij})^4}{r_{ij}^8} - \frac{\alpha (x_{ij}y_{ij}z_{ij})^4}{r_{ij}^{12}} \right] + b \right\} }{ \sum_{i \ne j} \sigma(r_{ij}) }
\f]

In this expression \f$x_{ij}\f$, \f$y_{ij}\f$ and \f$z_{ij}\f$ are the \f$x\f$, \f$y\f$ and \f$z\f$ components of the vector connecting atom \f$i\f$ to
atom \f$j\f$ and \f$r_{ij}\f$ is the magnitude of this vector.  \f$\sigma(r_{ij})\f$ is a \ref switchingfunction that acts on the distance between
atom \f$i\f$ and atom \f$j\f$ and its inclusion in the numerator and the denominator of the above expression as well as the fact that we are summing
over all of the other atoms in the system ensures that we are calculating an average
of the function of \f$x_{ij}\f$, \f$y_{ij}\f$ and \f$z_{ij}\f$ for the atoms in the first coordination sphere around atom \f$i\f$.  Lastly, \f$\alpha\f$
is a parameter that can be set by the user, which by default is equal to three.  The values of \f$a\f$ and \f$b\f$ are calculated from \f$\alpha\f$ using:

\f[
a = \frac{ 80080}{ 2717 + 16 \alpha} \qquad \textrm{and} \qquad b = \frac{ 16(\alpha - 143) }{2717 + 16\alpha}
\f]

This quantity is once again a multicolvar so you can compute it for multiple atoms using a single PLUMED action and then compute
the average value for the atoms in your system, the number of atoms that have an \f$s_i\f$ value that is more that some target and
so on.  Notice also that you can rotate the reference frame if you are using a non-standard unit cell.

\par Examples

The following input calculates the FCCUBIC parameter for the 64 atoms in the system
and then calculates and prints the average value for this quantity.

\plumedfile
FCCUBIC SPECIES=1-64 SWITCH={RATIONAL D_0=3.0 R_0=1.5} MEAN LABEL=d
PRINT ARG=d.* FILE=colv
\endplumedfile

*/
//+ENDPLUMEDOC


class Fccubic : public CubicHarmonicBase {
private:
  double alpha, a1, b1;
public:
  static void registerKeywords( Keywords& keys );
  explicit Fccubic(const ActionOptions&);
  double calculateCubicHarmonic( const Vector& distance, const double& d2, Vector& myder ) const override;
};

PLUMED_REGISTER_ACTION(Fccubic,"FCCUBIC")

void Fccubic::registerKeywords( Keywords& keys ) {
  CubicHarmonicBase::registerKeywords( keys );
  keys.add("compulsory","ALPHA","3.0","The alpha parameter of the angular function");
}

Fccubic::Fccubic(const ActionOptions&ao):
  Action(ao),
  CubicHarmonicBase(ao)
{
  // Scaling factors such that '1' corresponds to fcc lattice
  // and '0' corresponds to isotropic (liquid)
  parse("ALPHA",alpha);
  a1 = 80080. / (2717. + 16*alpha); b1 = 16.*(alpha-143)/(2717+16*alpha);
  log.printf("  setting alpha parameter equal to %f \n",alpha);
  // And setup the ActionWithVessel
  checkRead();
}

double Fccubic::calculateCubicHarmonic( const Vector& distance, const double& d2, Vector& myder ) const {
  double x2 = distance[0]*distance[0];
  double x4 = x2*x2;

  double y2 = distance[1]*distance[1];
  double y4 = y2*y2;

  double z2 = distance[2]*distance[2];
  double z4 = z2*z2;

  double r8 = pow( d2, 4 );
  double r12 = pow( d2, 6 );

  double tmp = ((x4*y4)+(x4*z4)+(y4*z4))/r8-alpha*x4*y4*z4/r12;

  double t0 = (x2*y4+x2*z4)/r8-alpha*x2*y4*z4/r12;
  double t1 = (y2*x4+y2*z4)/r8-alpha*y2*x4*z4/r12;
  double t2 = (z2*x4+z2*y4)/r8-alpha*z2*x4*y4/r12;
  double t3 = (2*tmp-alpha*x4*y4*z4/r12)/d2;

  myder[0]=4*a1*distance[0]*(t0-t3);
  myder[1]=4*a1*distance[1]*(t1-t3);
  myder[2]=4*a1*distance[2]*(t2-t3);

  return a1*tmp+b1;
}

}
}

