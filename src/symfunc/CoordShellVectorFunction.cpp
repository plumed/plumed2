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
#include "SymmetryFunctionBase.h"
#include "tools/LeptonCall.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"

#include <string>
#include <cmath>

using namespace std;

namespace PLMD {
namespace symfunc {

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

//+PLUMEDOC MCOLVAR TETRAHEDRAL
/*
Calculate the degree to which the environment about ions has a tetrahedral order.

We can measure the degree to which the atoms in the first coordination shell around any atom, \f$i\f$ is
is arranged like a tetrahedron using the following function.

\f[
 s(i) = \frac{1}{\sum_j \sigma( r_{ij} )} \sum_j \sigma( r_{ij} )\left[ \frac{(x_{ij} + y_{ij} + z_{ij})^3}{r_{ij}^3} +
                                                                        \frac{(x_{ij} - y_{ij} - z_{ij})^3}{r_{ij}^3} +
                                                                        \frac{(-x_{ij} + y_{ij} - z_{ij})^3}{r_{ij}^3} +
                                                                        \frac{(-x_{ij} - y_{ij} + z_{ij})^3}{r_{ij}^3} \right]
\f]

Here \f$r_{ij}\f$ is the magnitude of the vector connecting atom \f$i\f$ to atom \f$j\f$ and \f$x_{ij}\f$, \f$y_{ij}\f$ and \f$z_{ij}\f$
are its three components.  The function  \f$\sigma( r_{ij} )\f$ is a \ref switchingfunction that acts on the distance between
atoms \f$i\f$ and \f$j\f$.  The parameters of this function should be set so that the function is equal to one
when atom \f$j\f$ is in the first coordination sphere of atom \f$i\f$ and is zero otherwise.

\par Examples

The following command calculates the average value of the TETRAHEDRAL parameter for a set of 64 atoms all of the same type
and outputs this quantity to a file called colvar.

\plumedfile
tt: TETRAHEDRAL SPECIES=1-64 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MEAN
PRINT ARG=tt.mean FILE=colvar
\endplumedfile

The following command calculates the number of TETRAHEDRAL parameters that are greater than 0.8 in a set of 10 atoms.
In this calculation it is assumed that there are two atom types A and B and that the first coordination sphere of the
10 atoms of type A contains atoms of type B.  The formula above is thus calculated for ten different A atoms and within
it the sum over \f$j\f$ runs over 40 atoms of type B that could be in the first coordination sphere.

\plumedfile
tt: TETRAHEDRAL SPECIESA=1-10 SPECIESB=11-40 SWITCH={RATIONAL D_0=1.3 R_0=0.2} MORE_THAN={RATIONAL R_0=0.8}
PRINT ARG=tt.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC


class CoordShellVectorFunction : public SymmetryFunctionBase {
private: 
  LeptonCall function; 
public:
  static void registerKeywords( Keywords& keys );
  explicit CoordShellVectorFunction(const ActionOptions&);
  void compute( const double& val, const Vector& dir, MultiValue& myvals ) const override;
};

PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"COORDINATION_SHELL_FUNCTION")

void CoordShellVectorFunction::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::registerKeywords( keys );
  keys.add("compulsory","FUNCTION","the function of the bond vectors you would like to compute");
}

CoordShellVectorFunction::CoordShellVectorFunction(const ActionOptions&ao):
  Action(ao),
  SymmetryFunctionBase(ao)
{
  std::string myfunc; parse("FUNCTION",myfunc); std::vector<std::string> argname(4);
  argname[0]="x"; argname[1]="y"; argname[2]="z"; argname[3]="r"; function.set( myfunc, argname, this );
  addValueWithDerivatives(); checkRead();
}

void CoordShellVectorFunction::compute( const double& val, const Vector& distance, MultiValue& myvals ) const {
  std::vector<double> args(4); Vector myder; args[0]=distance[0]; args[1]=distance[1]; args[2]=distance[2]; args[3]=distance.modulo();
  double tmp = function.evaluate( args ); addToValue( 0, val*tmp, myvals ); double rder = function.evaluateDeriv( 3, args );
  for(unsigned i=0;i<3;++i) myder[i] = function.evaluateDeriv( i, args ) + rder*distance[i] / args[3];
  addWeightDerivative( 0, tmp, myvals ); addVectorDerivatives( 0, val*myder, myvals );
}

class CoordShellVectorFunctionShortcut : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit CoordShellVectorFunctionShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CoordShellVectorFunctionShortcut,"TETRAHEDRAL")
PLUMED_REGISTER_ACTION(CoordShellVectorFunctionShortcut,"SIMPLECUBIC")
PLUMED_REGISTER_ACTION(CoordShellVectorFunctionShortcut,"COORDINATION_SHELL_AVERAGE")

void CoordShellVectorFunctionShortcut::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::shortcutKeywords( keys );
  keys.add("compulsory","FUNCTION","the function of the bond vectors that you would like to evaluate"); 
} 

CoordShellVectorFunctionShortcut::CoordShellVectorFunctionShortcut(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::string myfunc;
  if( getName()=="TETRAHEDRAL" ) myfunc = "((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3"; 
  else if( getName()=="SIMPLECUBIC" ) myfunc = "(x^4+y^4+z^4)/(r^4)";
  else parse("FUNCTION",myfunc);
  std::string sp_str, specA, specB; parse("SPECIES",sp_str); parse("SPECIESA",specA); parse("SPECIESB",specB);
  std::map<std::string,std::string> keymap; multicolvar::MultiColvarBase::readShortcutKeywords( keymap, this ); 
  if( sp_str.length()==0 && specA.length()==0 ) error("no atoms were specified in input"); 
  SymmetryFunctionBase::expandMatrix( true, getShortcutLabel(),  sp_str, specA, specB, this );
  readInputLine( getShortcutLabel() + ": COORDINATION_SHELL_FUNCTION FUNCTION=" + myfunc + " WEIGHT=" + getShortcutLabel() + "_mat.w "
                    + "VECTORS1=" + getShortcutLabel() + "_mat.x VECTORS2=" + getShortcutLabel() + "_mat.y VECTORS3=" + getShortcutLabel() + "_mat.z " +
                    convertInputLineToString() );
  readInputLine( getShortcutLabel() + "_denom: COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_mat.w");
  readInputLine( getShortcutLabel() + "_n: MATHEVAL ARG1=" + getShortcutLabel() + " ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_n", "", keymap, this );
} 

}
}

