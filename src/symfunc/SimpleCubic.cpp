/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "SymmetryFunctionBase.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR SIMPLECUBIC
/*
Calculate whether or not the coordination spheres of atoms are arranged as they would be in a simple
cubic structure.

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


class SimpleCubic : public SymmetryFunctionBase {
public:
  static void shortcutKeywords( Keywords& keys );
  static void expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                              const std::map<std::string,std::string>& keys,
                              std::vector<std::vector<std::string> >& actions );
  static void registerKeywords( Keywords& keys );
  explicit SimpleCubic(const ActionOptions&);
  void compute( const double& val, const Vector& dir, MultiValue& myvals ) const ;
};

PLUMED_REGISTER_ACTION(SimpleCubic,"SIMPLECUBIC")
PLUMED_REGISTER_SHORTCUT(SimpleCubic,"SIMPLECUBIC")

void SimpleCubic::shortcutKeywords( Keywords& keys ) {
  SymmetryFunctionBase::shortcutKeywords( keys );
}

void SimpleCubic::expandShortcut( const std::string& lab, const std::vector<std::string>& words,
                                  const std::map<std::string,std::string>& keys,
                                  std::vector<std::vector<std::string> >& actions ) {
  SymmetryFunctionBase::expandMatrix( true, lab, words, keys, actions );
  // Input for simple cubic action
  std::vector<std::string> input; input.push_back(lab + ":"); input.push_back("SIMPLECUBIC");
  input.push_back("WEIGHT=" + lab + "_mat.w" ); input.push_back("VECTORS1=" + lab + "_mat.x" ); 
  input.push_back("VECTORS2=" + lab + "_mat.y" ); input.push_back("VECTORS3=" + lab + "_mat.z" ); 
  actions.push_back( input );
  // Input for denominator (coord)
  std::vector<std::string> d_input; d_input.push_back(lab + "_denom:"); d_input.push_back("COORDINATIONNUMBER");
  d_input.push_back("WEIGHT=" + lab + "_mat.w"); actions.push_back( d_input );
  // Input for matheval action
  std::vector<std::string> me_input; me_input.push_back(lab + "_n:"); me_input.push_back("MATHEVAL"); 
  me_input.push_back("ARG1=" + lab); me_input.push_back("ARG2=" + lab + "_denom"); me_input.push_back("FUNC=x/y");
  me_input.push_back("PERIODIC=NO"); actions.push_back(me_input);
  multicolvar::MultiColvarBase::expandFunctions( lab, lab + "_n", "", words, keys, actions );
}

void SimpleCubic::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys );
}

SimpleCubic::SimpleCubic(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  addValueWithDerivatives(); checkRead();
}

void SimpleCubic::compute( const double& val, const Vector& dir, MultiValue& myvals ) const {
  double x2 = dir[0]*dir[0];
  double x3 = dir[0]*x2;
  double x4 = dir[0]*x3;

  double y2 = dir[1]*dir[1];
  double y3 = dir[1]*y2;
  double y4 = dir[1]*y3;

  double z2 = dir[2]*dir[2];
  double z3 = dir[2]*z2;
  double z4 = dir[2]*z3;

  double d2 = dir.modulo2();
  double r4 = pow( d2, 2 );
  double tmp = ( x4 + y4 + z4 ) / r4;

  double t1=(x2+y2+z2), t2=t1*t1, t3=(x4+y4+z4)/(t1*t2);

  Vector myder; 
  myder[0] = 4*x3/t2-4*dir[0]*t3;
  myder[1] = 4*y3/t2-4*dir[1]*t3;
  myder[2] = 4*z3/t2-4*dir[2]*t3;
  
  // Accumulate derivatives
  addToValue( 0, val*tmp, myvals );
  addWeightDerivative( 0, tmp, myvals ); 
  addVectorDerivatives( 0, val*myder, myvals ); 
}

}
}

