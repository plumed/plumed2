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
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include "core/ActionWithValue.h"
#include "multicolvar/MultiColvarShortcuts.h"

//+PLUMEDOC MCOLVAR ATOMIC_SMAC
/*
Calculate the atomic smac CV

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace symfunc {

class AtomicSMAC : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit AtomicSMAC(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(AtomicSMAC,"ATOMIC_SMAC")

void AtomicSMAC::registerKeywords(Keywords& keys) {
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
  keys.needsAction("CONTACT_MATRIX");
  keys.needsAction("GSYMFUNC_THREEBODY");
  keys.needsAction("ONES");
  keys.needsAction("MATRIX_VECTOR_PRODUCT");
}

AtomicSMAC::AtomicSMAC(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Create the matrices
  std::string sw_input;
  parse("SWITCH",sw_input);
  std::string sp_lab, sp_laba;
  parse("SPECIES",sp_lab);
  parse("SPECIESA",sp_laba);
  std::string cmap_input = getShortcutLabel() + "_cmap: CONTACT_MATRIX";
  if( sp_lab.length()>0 ) {
    readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUP=" + sp_lab + " COMPONENTS SWITCH={" + sw_input + "}");
  } else if( sp_laba.length()>0 ) {
    std::string sp_labb;
    parse("SPECIESB",sp_labb);
    readInputLine( getShortcutLabel() + "_cmap: CONTACT_MATRIX GROUPA=" + sp_laba + " GROUPB=" + sp_labb + " COMPONENTS SWITCH={" + sw_input + "}");
  }
  // Now need the Gaussians
  std::string mykernels;
  for(unsigned i=1;; ++i) {
    std::string kstr_inpt, istr, kern_str;
    Tools::convert( i, istr );
    if( !parseNumbered("KERNEL",i,kstr_inpt ) ) {
      break;
    }
    std::vector<std::string> words = Tools::getWords(kstr_inpt);
    if( words[0]=="GAUSSIAN" ) {
      kern_str="gaussian";
    } else {
      error("unknown kernel type");
    }
    std::string center, var;
    Tools::parse(words,"CENTER",center);
    Tools::parse(words,"SIGMA",var);
    if( mykernels.length()==0 ) {
      mykernels = "exp(-(ajik-" + center + ")^2/(2*" + var + "*" + var + "))";
    } else {
      mykernels = mykernels + "+exp(-(ajik-" + center + ")^2/(2*" + var + "*" + var + "))";
    }
  }
  // Hard coded switching function on minimum distance here -- should be improved
  readInputLine( getShortcutLabel() + "_ksum: GSYMFUNC_THREEBODY WEIGHT=" + getShortcutLabel() + "_cmap.w " +
                 "ARG=" + getShortcutLabel() + "_cmap.x," + getShortcutLabel() + "_cmap.y," + getShortcutLabel() + "_cmap.z"
                 " FUNCTION1={FUNC=" + mykernels + " LABEL=n} FUNCTION2={FUNC=1 LABEL=d}" );
  // And just the sum of the coordination numbers
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_cmap");
  plumed_assert( av && av->getNumberOfComponents()>0 && (av->copyOutput(0))->getRank()==2 );
  std::string size;
  Tools::convert( (av->copyOutput(0))->getShape()[1], size );
  readInputLine( getShortcutLabel() + "_ones: ONES SIZE=" + size );
  readInputLine( getShortcutLabel() + "_denom: MATRIX_VECTOR_PRODUCT ARG=" + getShortcutLabel() + "_cmap.w," + getShortcutLabel() + "_ones");
  // And the transformed switching functions
  std::string swcoord_str;
  parse("SWITCH_COORD",swcoord_str);
  readInputLine( getShortcutLabel() + "_mtdenom: MORE_THAN ARG=" + getShortcutLabel() + "_denom SWITCH={" + swcoord_str +"}");
  // And matheval to get the final quantity
  readInputLine( getShortcutLabel() + "_smac: CUSTOM ARG=" + getShortcutLabel() + "_ksum.n," + getShortcutLabel() + "_mtdenom," + getShortcutLabel() + "_ksum.d FUNC=x*y/z PERIODIC=NO");
  // And this expands everything
  multicolvar::MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + "_smac", "", this );
}

}
}
