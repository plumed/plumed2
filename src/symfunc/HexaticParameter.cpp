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
#include "CoordinationNumbers.h"
#include "core/ActionShortcut.h"
#include "multicolvar/MultiColvarBase.h"
#include "core/ActionRegister.h"

#include <complex>

using namespace std;

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR HEXACTIC_PARAMETER
/*

\bug Virial is not working currently

\par Examples


*/
//+ENDPLUMEDOC


class HexacticParameter : public ActionShortcut {
private:
  void createVectorNormInput( const std::string& ilab, const std::string& olab, const std::string& vlab );
public:
  static void registerKeywords( Keywords& keys );
  explicit HexacticParameter(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(HexacticParameter,"HEXACTIC_PARAMETER")

void HexacticParameter::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.add("compulsory","PLANE","the plane to use when calculating the value of the order parameter should be xy, xz or yz");
  keys.addFlag("VMEAN",false,"calculate the norm of the mean vector.");
  keys.addOutputComponent("_vmean","VMEAN","the norm of the mean vector");
  keys.addFlag("VSUM",false,"calculate the norm of the sum of all the vectors");
  keys.addOutputComponent("_vsum","VSUM","the norm of the mean vector");
}

HexacticParameter::HexacticParameter( const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::string sp_str, specA, specB; parse("SPECIES",sp_str); parse("SPECIESA",specA); parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( true, getShortcutLabel(), sp_str, specA, specB, this ); 
  std::string myplane; parse("PLANE",myplane);
  if( myplane=="xy" ) {
     readInputLine( getShortcutLabel() + ": CYLINDRICAL_HARMONIC_MATRIX DEGREE=6 ARG1=" + getShortcutLabel() + "_mat.x ARG2=" + getShortcutLabel() + "_mat.y ARG3=" + getShortcutLabel() + "_mat.w" );
  } else if( myplane=="xz" ) {
     readInputLine( getShortcutLabel() + ": CYLINDRICAL_HARMONIC_MATRIX DEGREE=6 ARG1=" + getShortcutLabel() + "_mat.x ARG2=" + getShortcutLabel() + "_mat.z ARG3=" + getShortcutLabel() + "_mat.w" );
  } else if( myplane=="yz" ) {
     readInputLine( getShortcutLabel() + ": CYLINDRICAL_HARMONIC_MATRIX DEGREE=6 ARG1=" + getShortcutLabel() + "_mat.y ARG2=" + getShortcutLabel() + "_mat.z ARG3=" + getShortcutLabel() + "_mat.w" );
  } else error("invalid input for plane -- should be xy, xz or yz");
  // And coordination number  
  readInputLine( getShortcutLabel() + "_rm: COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + ".rm");
  readInputLine( getShortcutLabel() + "_im: DOT ARG1=" + getShortcutLabel() + ".im ARG2=" + getShortcutLabel() + "_rm_ones");
  // Input for denominator (coord)
  readInputLine( getShortcutLabel() + "_denom: DOT ARG1=" + getShortcutLabel() + "_mat.w ARG2=" + getShortcutLabel() + "_rm_ones");
  // Divide real part by coordination numbers
  readInputLine( getShortcutLabel() + "_rmn: CUSTOM ARG1=" + getShortcutLabel() + "_rm ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  // Devide imaginary part by coordination number
  readInputLine( getShortcutLabel() + "_imn: CUSTOM ARG1=" + getShortcutLabel() + "_im ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");

  // If we are doing VMEAN determine sum of vector components
  bool do_vmean; parseFlag("VMEAN",do_vmean);
  if( do_vmean ) {
    // Real part
    readInputLine( getShortcutLabel() + "_rms: MEAN ARG=" + getShortcutLabel() + "_rmn PERIODIC=NO");
    // Imaginary part
    readInputLine( getShortcutLabel() + "_ims: MEAN ARG=" + getShortcutLabel() + "_imn PERIODIC=NO");
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vmean", "ms" );
  }
  bool do_vsum; parseFlag("VSUM",do_vsum);
  if( do_vsum ) {
    // Real part
    readInputLine( getShortcutLabel() + "_rmz: SUM ARG=" + getShortcutLabel() + "_rmn PERIODIC=NO"); 
    // Imaginary part
    readInputLine( getShortcutLabel() + "_imz: SUM ARG=" + getShortcutLabel() + "_imn PERIODIC=NO"); 
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vsum", "mz" );
  }

  // Now calculate the total length of the vector
  createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_norm", "mn" );
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_norm", "", this );
}

void HexacticParameter::createVectorNormInput( const std::string& ilab, const std::string& olab, const std::string& vlab ) {
  readInputLine( olab + "2: COMBINE PERIODIC=NO ARG1=" + ilab + "_r" + vlab + " ARG2=" + ilab + "_i" + vlab + " POWERS=2,2" );
  readInputLine( olab + ": MATHEVAL ARG1=" + olab + "2 FUNC=sqrt(x) PERIODIC=NO"); 
}

}
}

