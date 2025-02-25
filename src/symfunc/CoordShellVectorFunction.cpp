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
#include "CoordinationNumbers.h"
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/ActionWithValue.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR COORDINATION_SHELL_FUNCTION
/*
Calculate an arbitrary function of all the bond vectors in the first coordination sphere of an atom

\par Examples


*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR COORDINATION_SHELL_AVERAGE
/*
Calculate an arbitrary function of all the bond vectors in the first coordination sphere of an atom and take an average

\par Examples


*/
//+ENDPLUMEDOC

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

class CoordShellVectorFunction : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit CoordShellVectorFunction(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"FCCUBIC")
PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"TETRAHEDRAL")
PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"SIMPLECUBIC")
PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"COORDINATION_SHELL_FUNCTION")
PLUMED_REGISTER_ACTION(CoordShellVectorFunction,"COORDINATION_SHELL_AVERAGE")

void CoordShellVectorFunction::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.add("compulsory","FUNCTION","the function of the bond vectors that you would like to evaluate");
  keys.add("compulsory","PHI","0.0","The Euler rotational angle phi");
  keys.add("compulsory","THETA","0.0","The Euler rotational angle theta");
  keys.add("compulsory","PSI","0.0","The Euler rotational angle psi");
  keys.add("compulsory","ALPHA","3.0","The alpha parameter of the angular function that is used for FCCUBIC");
  keys.addFlag("LOWMEM",false,"this flag does nothing and is present only to ensure back-compatibility");
  keys.setValueDescription("vector","the symmetry function for each of the specified atoms");
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("FCCUBIC_FUNC");
  keys.needsAction("CUSTOM");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
}

CoordShellVectorFunction::CoordShellVectorFunction(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string matlab, sp_str, specA, specB;
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
  parse("SPECIES",sp_str);
  parse("SPECIESA",specA);
  parse("SPECIESB",specB);
  if( sp_str.length()>0 || specA.length()>0 ) {
    matlab = getShortcutLabel() + "_mat";
    CoordinationNumbers::expandMatrix( true, getShortcutLabel(),  sp_str, specA, specB, this );
  } else {
    error("found no input atoms use SPECIES/SPECIESA");
  }
  double phi, theta, psi;
  parse("PHI",phi);
  parse("THETA",theta);
  parse("PSI",psi);
  std::vector<std::string> rotelements(9);
  std::string xvec = matlab + ".x", yvec = matlab + ".y", zvec = matlab + ".z";
  if( phi!=0 || theta!=0 || psi!=0 ) {
    Tools::convert( std::cos(psi)*std::cos(phi)-std::cos(theta)*std::sin(phi)*std::sin(psi), rotelements[0] );
    Tools::convert( std::cos(psi)*std::sin(phi)+std::cos(theta)*std::cos(phi)*std::sin(psi), rotelements[1] );
    Tools::convert( std::sin(psi)*std::sin(theta), rotelements[2] );

    Tools::convert( -std::sin(psi)*std::cos(phi)-std::cos(theta)*std::sin(phi)*std::cos(psi), rotelements[3] );
    Tools::convert( -std::sin(psi)*std::sin(phi)+std::cos(theta)*std::cos(phi)*std::cos(psi), rotelements[4] );
    Tools::convert( std::cos(psi)*std::sin(theta), rotelements[5] );

    Tools::convert( std::sin(theta)*std::sin(phi), rotelements[6] );
    Tools::convert( -std::sin(theta)*std::cos(phi), rotelements[7] );
    Tools::convert( std::cos(theta), rotelements[8] );
    readInputLine( getShortcutLabel() + "_xrot: CUSTOM ARG=" + matlab + ".x," + matlab + ".y," + matlab + ".z FUNC=" + rotelements[0] + "*x+" + rotelements[1] + "*y+" + rotelements[2] + "*z PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_yrot: CUSTOM ARG=" + matlab + ".x," + matlab + ".y," + matlab + ".z FUNC=" + rotelements[3] + "*x+" + rotelements[4] + "*y+" + rotelements[5] + "*z PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_zrot: CUSTOM ARG=" + matlab + ".x," + matlab + ".y," + matlab + ".z FUNC=" + rotelements[6] + "*x+" + rotelements[7] + "*y+" + rotelements[8] + "*z PERIODIC=NO");
  }
  // Calculate FCC cubic function from bond vectors
  if( getName()=="FCCUBIC" ) {
    std::string alpha;
    parse("ALPHA",alpha);
    readInputLine( getShortcutLabel() + "_vfunc: FCCUBIC_FUNC ARG=" + xvec + "," + yvec + "," + zvec+ " ALPHA=" + alpha);
  } else if( getName()=="TETRAHEDRAL" ) {
    readInputLine( getShortcutLabel() + "_r: CUSTOM ARG=" + xvec + "," + yvec + "," + zvec + " PERIODIC=NO FUNC=sqrt(x*x+y*y+z*z)");
    readInputLine( getShortcutLabel() + "_vfunc: CUSTOM ARG=" + xvec + "," + yvec + "," + zvec + "," + getShortcutLabel() + "_r"
                   + " VAR=x,y,z,r PERIODIC=NO FUNC=((x+y+z)/r)^3+((x-y-z)/r)^3+((-x+y-z)/r)^3+((-x-y+z)/r)^3" );
  } else if( getName()=="SIMPLECUBIC" ) {
    readInputLine( getShortcutLabel() + "_r: CUSTOM ARG=" + xvec + "," + yvec + "," + zvec + " PERIODIC=NO FUNC=sqrt(x*x+y*y+z*z)");
    readInputLine( getShortcutLabel() + "_vfunc: CUSTOM ARG=" + xvec + "," + yvec + "," + zvec + "," + getShortcutLabel() + "_r"
                   + " VAR=x,y,z,r PERIODIC=NO FUNC=(x^4+y^4+z^4)/(r^4)" );
  } else {
    std::string myfunc;
    parse("FUNCTION",myfunc);
    if( myfunc.find("r")!=std::string::npos ) {
      readInputLine( getShortcutLabel() + "_r: CUSTOM ARG=" + xvec + "," + yvec + "," + zvec + " PERIODIC=NO FUNC=sqrt(x*x+y*y+z*z)");
      readInputLine( getShortcutLabel() + "_vfunc: CUSTOM ARG=" + xvec + "," + yvec + "," + zvec + "," + getShortcutLabel() + "_r VAR=x,y,z,r PERIODIC=NO FUNC=" + myfunc );
    } else {
      readInputLine( getShortcutLabel() + "_vfunc: CUSTOM ARG=" + xvec + "," + yvec + "," + zvec + " PERIODIC=NO FUNC=" + myfunc );
    }
  }
  // Hadamard product of function above and weights
  readInputLine( getShortcutLabel() + "_wvfunc: CUSTOM ARG=" + getShortcutLabel() + "_vfunc," + matlab + ".w FUNC=x*y PERIODIC=NO");
  // And coordination numbers
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_wvfunc," + getShortcutLabel() + "_ones");
  std::string olab=getShortcutLabel();
  if( getName()!="COORDINATION_SHELL_FUNCTION" ) {
    olab = getShortcutLabel() + "_n";
    // Calculate coordination numbers for denominator
    readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + matlab + ".w," + getShortcutLabel() + "_ones");
    // And normalise
    readInputLine( getShortcutLabel() + "_n: CUSTOM ARG=" + getShortcutLabel() + "," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  }
  // And expand the functions
  std::map<std::string,std::string> keymap;
  multicolvar::MultiColvarShortcuts::readShortcutKeywords( keymap, this );
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), olab, "", keymap, this );
}

}
}

