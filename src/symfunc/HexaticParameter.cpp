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
#include "multicolvar/MultiColvarShortcuts.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionRegister.h"
#include "core/ActionWithValue.h"
#include "CoordinationNumbers.h"

#include <complex>

namespace PLMD {
namespace symfunc {

//+PLUMEDOC MCOLVAR HEXACTIC_PARAMETER
/*
Calculate the hexatic order parameter

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
  keys.needsAction("CYLINDRICAL_HARMONIC_MATRIX"); keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT"); keys.needsAction("CUSTOM");
  keys.needsAction("MEAN"); keys.needsAction("SUM"); keys.needsAction("COMBINE");
}

HexacticParameter::HexacticParameter( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string sp_str, specA, specB; parse("SPECIES",sp_str); parse("SPECIESA",specA); parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( true, getShortcutLabel(), sp_str, specA, specB, this );
  std::string myplane; parse("PLANE",myplane);
  if( myplane=="xy" ) {
    readInputLine( getShortcutLabel() + ": CYLINDRICAL_HARMONIC_MATRIX DEGREE=6 ARG=" + getShortcutLabel() + "_mat.x," + getShortcutLabel() + "_mat.y," + getShortcutLabel() + "_mat.w" );
  } else if( myplane=="xz" ) {
    readInputLine( getShortcutLabel() + ": CYLINDRICAL_HARMONIC_MATRIX DEGREE=6 ARG=" + getShortcutLabel() + "_mat.x," + getShortcutLabel() + "_mat.z," + getShortcutLabel() + "_mat.w" );
  } else if( myplane=="yz" ) {
    readInputLine( getShortcutLabel() + ": CYLINDRICAL_HARMONIC_MATRIX DEGREE=6 ARG=" + getShortcutLabel() + "_mat.y," + getShortcutLabel() + "_mat.z," + getShortcutLabel() + "_mat.w" );
  } else error("invalid input for plane -- should be xy, xz or yz");
  // And coordination number
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_mat");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size; Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_rm: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + ".rm," + getShortcutLabel() + "_ones");
  readInputLine( getShortcutLabel() + "_im: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + ".im," + getShortcutLabel() + "_ones");
  // Input for denominator (coord)
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_mat.w," + getShortcutLabel() + "_ones");
  // Divide real part by coordination numbers
  readInputLine( getShortcutLabel() + "_rmn: CUSTOM ARG=" + getShortcutLabel() + "_rm," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  // Devide imaginary part by coordination number
  readInputLine( getShortcutLabel() + "_imn: CUSTOM ARG=" + getShortcutLabel() + "_im," + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");

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
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_norm", "", this );
}

void HexacticParameter::createVectorNormInput( const std::string& ilab, const std::string& olab, const std::string& vlab ) {
  readInputLine( olab + "2: COMBINE PERIODIC=NO ARG=" + ilab + "_r" + vlab + "," + ilab + "_i" + vlab + " POWERS=2,2" );
  readInputLine( olab + ": CUSTOM ARG=" + olab + "2 FUNC=sqrt(x) PERIODIC=NO");
}

}
}

