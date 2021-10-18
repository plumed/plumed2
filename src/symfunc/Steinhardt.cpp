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

//+PLUMEDOC MCOLVAR SPHERICAL_HARMONIC
/*

\par Examples


*/
//+ENDPLUMEDOC


class Steinhardt : public ActionShortcut {
private: 
  void createVectorNormInput( const std::string& ilab, const std::string& olab, const int& l, const std::string& sep, const std::string& vlab );
public:
  static void registerKeywords( Keywords& keys );
  explicit Steinhardt(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Steinhardt,"Q1")
PLUMED_REGISTER_ACTION(Steinhardt,"Q3")
PLUMED_REGISTER_ACTION(Steinhardt,"Q4")
PLUMED_REGISTER_ACTION(Steinhardt,"Q6")

void Steinhardt::registerKeywords( Keywords& keys ) {
  CoordinationNumbers::shortcutKeywords( keys );
  keys.addFlag("VMEAN",false,"calculate the norm of the mean vector.");
  keys.addOutputComponent("_vmean","VMEAN","the norm of the mean vector");
  keys.addFlag("VSUM",false,"calculate the norm of the sum of all the vectors");
  keys.addOutputComponent("_vsum","VSUM","the norm of the mean vector");
}

Steinhardt::Steinhardt( const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  std::string sp_str, specA, specB; parse("SPECIES",sp_str); parse("SPECIESA",specA); parse("SPECIESB",specB);
  CoordinationNumbers::expandMatrix( true, getShortcutLabel(), sp_str, specA, specB, this ); int l;
  std::string sph_input = getShortcutLabel() + "_sh: SPHERICAL_HARMONIC ARG1=" + getShortcutLabel() + "_mat.x ARG2=" + getShortcutLabel() + "_mat.y ARG3=" + getShortcutLabel() + "_mat.z ARG4=" + getShortcutLabel() + "_mat.w";

  if( getName()=="Q1" ) {
    sph_input +=" L=1"; l=1;
  } else if( getName()=="Q3" ) {
    sph_input += " L=3"; l=3;
  } else if( getName()=="Q4" ) {
    sph_input += " L=4"; l=4;
  } else if( getName()=="Q6" ) {
    sph_input += " L=6"; l=6;
  } else {
    plumed_merror("invalid input");
  }
  readInputLine( sph_input );

  // Input for denominator (coord)
  readInputLine( getShortcutLabel() + "_denom: COORDINATIONNUMBER WEIGHT=" + getShortcutLabel() + "_mat.w");
  std::string snum; Tools::convert( -l, snum ); std::string arg2=getShortcutLabel() + "_denom_ones," + getShortcutLabel() + "_denom_ones"; 
  std::string arg1=getShortcutLabel() + "_sh.rm-[" + snum + "]," + getShortcutLabel() + "_sh.im-[" + snum + "]";
  for(int i=-l+1; i<=l; ++i) { 
      Tools::convert( i, snum ); arg1+= "," + getShortcutLabel() + "_sh.rm-[" + snum + "]," + getShortcutLabel() + "_sh.im-[" + snum + "]"; 
      arg2 += "," + getShortcutLabel() + "_denom_ones," + getShortcutLabel() + "_denom_ones";
  }
  readInputLine( getShortcutLabel() + ": DOT ARG1=" + arg1 + " ARG2=" + arg2 );

  // If we are doing VMEAN determine sum of vector components
  bool do_vmean; parseFlag("VMEAN",do_vmean);
  bool do_vsum; parseFlag("VSUM",do_vsum);
  if( do_vmean || do_vsum ) {
    // Divide all components by coordination numbers
    for(int i=-l; i<=l; ++i) {
      Tools::convert( i, snum );
      // Real part
      readInputLine( getShortcutLabel() + "_rmn-[" + snum + "]: CUSTOM ARG1=" + getShortcutLabel() + ".rm-[" + snum + "] ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
      // Imaginary part
      readInputLine( getShortcutLabel() + "_imn-[" + snum + "]: CUSTOM ARG1=" + getShortcutLabel() + ".im-[" + snum + "] ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
    }
  }

  if( do_vmean ) { 
    for(int i=-l; i<=l; ++i) {
      std::string snum; Tools::convert( i, snum );
      // Real part
      readInputLine( getShortcutLabel() + "_rms-[" + snum + "]: MEAN ARG=" + getShortcutLabel() + "_rmn-[" + snum + "] PERIODIC=NO");
      // Imaginary part
      readInputLine( getShortcutLabel() + "_ims-[" + snum + "]: MEAN ARG=" + getShortcutLabel() + "_imn-[" + snum + "] PERIODIC=NO");
    }
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vmean", l, "_", "ms" );
  }
  if( do_vsum ) {
    for(int i=-l; i<=l; ++i) {
      std::string snum; Tools::convert( i, snum );
      // Real part
      readInputLine( getShortcutLabel() + "_rmz-[" + snum + "]: SUM ARG=" + getShortcutLabel() + "_rmn-[" + snum + "] PERIODIC=NO"); 
      // Imaginary part
      readInputLine( getShortcutLabel() + "_imz-[" + snum + "]: SUM ARG=" + getShortcutLabel() + "_imn-[" + snum + "] PERIODIC=NO"); 
    }
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vsum", l, "_", "mz" );
  }

  // Now calculate the total length of the vector
  createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_norm", l, ".", "m" );
  // And take average
  readInputLine( getShortcutLabel() + "_anorm: CUSTOM ARG1=" + getShortcutLabel() + "_norm ARG2=" + getShortcutLabel() + "_denom FUNC=x/y PERIODIC=NO");
  multicolvar::MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_anorm", "", this );
}

void Steinhardt::createVectorNormInput( const std::string& ilab, const std::string& olab, const int& l, const std::string& sep, const std::string& vlab ) {
  std::string arg_inp, norm_input = olab + "2: COMBINE PERIODIC=NO POWERS=2"; arg_inp = "";
  std::string snum, num; unsigned nn=1;
  for(int i=-l; i<=l; ++i) {
    Tools::convert( nn, num ); Tools::convert( i, snum );
    arg_inp += " ARG" + num + "=" + ilab + sep + "r" + vlab + "-[" + snum + "]";
    nn++; Tools::convert( nn, num );
    arg_inp += " ARG" + num + "=" + ilab + sep + "i" + vlab + "-[" + snum + "]";
    nn++;
    if( i==-l ) norm_input += ",2"; else norm_input += ",2,2";
  }
  readInputLine( norm_input + arg_inp );
  readInputLine( olab + ": MATHEVAL ARG1=" + olab + "2 FUNC=sqrt(x) PERIODIC=NO"); 
}

}
}

