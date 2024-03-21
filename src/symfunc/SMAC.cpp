/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2023 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "multicolvar/MultiColvarShortcuts.h"

//+PLUMEDOC MCOLVAR SMAC
/*
Calculate the SMAC order parameter for a set of molecules

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace symfunc {

class SMAC : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit SMAC(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(SMAC,"SMAC")

void SMAC::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","SPECIES","");
  keys.add("optional","SPECIESA","");
  keys.add("optional","SPECIESB","");
  keys.add("optional","SWITCH","This keyword is used if you want to employ an alternative to the continuous swiching function defined above. "
           "The following provides information on the \\ref switchingfunction that are available. "
           "When this keyword is present you no longer need the NN, MM, D_0 and R_0 keywords.");
  keys.add("numbered","KERNEL","The kernels used in the function of the angle");
  keys.add("optional","SWITCH_COORD","This keyword is used to define the coordination switching function.");
  keys.reset_style("KERNEL","optional");
  multicolvar::MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("VSTACK"); keys.needsAction("TRANSPOSE"); keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("TORSIONS_MATRIX"); keys.needsAction("COMBINE"); keys.needsAction("CUSTOM");
  keys.needsAction("ONES"); keys.needsAction("MATRIX_VECTOR_PRODUCT"); keys.needsAction("MORE_THAN");
}

SMAC::SMAC(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  // Create the matrices
  std::string sw_input; parse("SWITCH",sw_input);
  std::string sp_lab, sp_laba; parse("SPECIES",sp_lab); parse("SPECIESA",sp_laba);
  if( sp_lab.length()>0 ) {
    readInputLine( getShortcutLabel() + "_vecs: VSTACK ARG=" + sp_lab + ".x," + sp_lab + ".y," + sp_lab + ".z" );
    readInputLine( getShortcutLabel() + "_vecsT: TRANSPOSE ARG=" + getShortcutLabel() + "_vecs" );
    readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUP=" + sp_lab + " SWITCH={" + sw_input + "}");
    readInputLine( getShortcutLabel() + "_tpmat: TORSIONS_MATRIX ARG=" + getShortcutLabel() + "_vecs," + getShortcutLabel() + "_vecsT POSITIONS1=" + sp_lab + " POSITIONS2=" + sp_lab );
  } else if( sp_laba.length()>0 ) {
    std::string sp_labb; parse("SPECIESB",sp_labb);
    readInputLine( getShortcutLabel() + "_vecsa: VSTACK ARG=" + sp_laba + ".x," + sp_laba + ".y," + sp_laba + ".z" );
    readInputLine( getShortcutLabel() + "_vecsb: VSTACK ARG=" + sp_labb + ".x," + sp_labb + ".y," + sp_labb + ".z" );
    readInputLine( getShortcutLabel() + "_vecsbT: TRANSPOSE ARG=" + getShortcutLabel() + "_vecsb" );
    readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUPA=" + sp_laba + " GROUPB=" + sp_labb + " SWITCH={" + sw_input + "}");
    readInputLine( getShortcutLabel() + "_tpmat: TORSIONS_MATRIX ARG=" + getShortcutLabel() + "_vecsa," + getShortcutLabel() + "_vecsbT POSITIONS1=" + sp_laba + " POSITIONS2=" + sp_labb );
  }
  // Now need the Gaussians
  std::string kmap_input= getShortcutLabel() + "_ksum: COMBINE PERIODIC=NO";
  for(unsigned i=1;; ++i) {
    std::string kstr_inpt, istr; Tools::convert( i, istr );
    if( !parseNumbered("KERNEL",i,kstr_inpt ) ) { break; }
    std::vector<std::string> words = Tools::getWords(kstr_inpt);
    std::string center, var; Tools::parse(words,"CENTER",center); Tools::parse(words,"SIGMA",var);
    double nsig; Tools::convert( var, nsig ); std::string coeff; Tools::convert( 1/(nsig*nsig), coeff );
    readInputLine( getShortcutLabel() + "_kf" + istr + "_r2: COMBINE PERIODIC=NO ARG=" + getShortcutLabel() + "_tpmat COEFFICIENTS=" + coeff + " PARAMETERS=" + center + " POWERS=2");
    if( words[0]=="GAUSSIAN" ) readInputLine( getShortcutLabel() + "_kf" + istr + ": CUSTOM PERIODIC=NO FUNC=exp(-x/2) ARG=" +  getShortcutLabel() + "_kf" + istr + "_r2" );
    else if( words[0]=="TRIANGULAR" ) readInputLine( getShortcutLabel() + "_kf" + istr + ": CUSTOM PERIODIC=NO FUNC=step(1-sqrt(x))*(1-sqrt(x)) ARG=" + getShortcutLabel() + "_kf" + istr + "_r2" );
    else readInputLine( getShortcutLabel() + "_kf" + istr + ": CUSTOM PERIODIC=NO FUNC=" + words[0] + " ARG=" + getShortcutLabel() + "_kf" + istr + "_r2" );
    if( i==1 ) kmap_input += " ARG=" + getShortcutLabel() + "_kf" + istr;
    else kmap_input += "," + getShortcutLabel() + "_kf" + istr;
  }
  readInputLine( kmap_input );
  // Now create the product matrix
  readInputLine( getShortcutLabel() + "_prod: CUSTOM ARG=" + getShortcutLabel() + "_cmap," + getShortcutLabel() + "_ksum FUNC=x*y PERIODIC=NO");
  // Now the sum of coordination numbers times the switching functions
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_cmap");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size; Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + ": MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_prod," + getShortcutLabel() + "_ones");
  // And just the sum of the coordination numbers
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_cmap," + getShortcutLabel() + "_ones");
  // And the transformed switching functions
  std::string swcoord_str; parse("SWITCH_COORD",swcoord_str);
  readInputLine( getShortcutLabel() + "_mtdenom: MORE_THAN ARG=" + getShortcutLabel() + "_denom SWITCH={" + swcoord_str +"}");
// And matheval to get the final quantity
  readInputLine( getShortcutLabel() + "_smac: CUSTOM ARG=" + getShortcutLabel() + "," + getShortcutLabel() + "_mtdenom," + getShortcutLabel() + "_denom FUNC=(x*y)/z PERIODIC=NO");
  // And this expands everything
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_smac", "", this );
}

}
}
