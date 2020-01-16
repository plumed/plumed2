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
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "SymmetryFunctionBase.h"
#include "multicolvar/MultiColvarBase.h"

#include <complex>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR LOCAL_CRYSTALINITY
/*

\par Examples


*/
//+ENDPLUMEDOC


class LocalCrystallinity : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit LocalCrystallinity(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(LocalCrystallinity,"LOCAL_CRYSTALINITY")

void LocalCrystallinity::registerKeywords( Keywords& keys ) {
  SymmetryFunctionBase::shortcutKeywords( keys );
  keys.add("numbered","GVECTOR","the coefficients of the linear combinations to compute for the CV");
}

LocalCrystallinity::LocalCrystallinity( const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // This builds an adjacency matrix
  std::string sp_str, specA, specB; parse("SPECIES",sp_str); parse("SPECIESA",specA); parse("SPECIESB",specB);
  SymmetryFunctionBase::expandMatrix( true, getShortcutLabel(), sp_str, specA, specB, this ); 
  // Input for denominator (coord)
  readInputLine( getShortcutLabel() + "_denom: COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_mat.w");
  // Input for numerator
  std::string finput = getShortcutLabel() + ": COMBINE NORMALIZE PERIODIC=NO";
  for(unsigned i=1;; ++i) {
    std::vector<std::string> gvec; std::string istr; Tools::convert( i, istr );
    if( !parseNumberedVector("GVECTOR",i,gvec) ) { break; }
    if( gvec.size()!=3 ) error("gvectors should have size 3");
    // This is the dot product between the input gvector and the bond vector
    readInputLine( getShortcutLabel() + "_dot" + istr + ": COMBINE ARG1=" + getShortcutLabel() + "_mat.x ARG2=" + getShortcutLabel() + 
                                                        "_mat.y ARG3=" + getShortcutLabel() + "_mat.z PERIODIC=NO COEFFICIENTS=" + gvec[0] + "," + gvec[1] + "," + gvec[2] );  
    // Now calculate the sine and cosine of the dot product
    readInputLine( getShortcutLabel() + "_cos" + istr + ": MATHEVAL ARG1=" + getShortcutLabel() +"_mat.w ARG2=" + getShortcutLabel() + "_dot" + istr + " FUNC=x*cos(y) PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_sin" + istr + ": MATHEVAL ARG1=" + getShortcutLabel() +"_mat.w ARG2=" + getShortcutLabel() + "_dot" + istr + " FUNC=x*sin(y) PERIODIC=NO");
    // And sum up the sine and cosine over the coordination spheres
    readInputLine( getShortcutLabel() + "_cossum" + istr + ": COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_cos" + istr );
    readInputLine( getShortcutLabel() + "_sinsum" + istr + ": COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_sin" + istr );
    // And average the sine and cosine over the number of bonds
    readInputLine( getShortcutLabel() + "_cosmean" + istr + ": MATHEVAL FUNC=x/y PERIODIC=NO ARG1=" + getShortcutLabel() + "_cossum" + istr + 
                                                                                           " ARG2=" + getShortcutLabel() + "_denom");
    readInputLine( getShortcutLabel() + "_sinmean" + istr + ": MATHEVAL FUNC=x/y PERIODIC=NO ARG1=" + getShortcutLabel() + "_sinsum" + istr + 
                                                                                           " ARG2=" + getShortcutLabel() + "_denom");
    // And work out the square modulus of this complex number
    readInputLine( getShortcutLabel() + "_" + istr + ": MATHEVAL FUNC=x*x+y*y PERIODIC=NO ARG1=" + getShortcutLabel() + "_cosmean" + istr + 
                                                                                        " ARG2=" + getShortcutLabel() + "_sinmean" + istr);
    // These are all the kvectors that we are adding together in the final combine for the final CV
    finput += " ARG" + istr + "=" + getShortcutLabel() + "_" + istr; 
  }
  // This computes the final CV
  readInputLine( finput );
  // Now calculate the total length of the vector
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
}

}
}

