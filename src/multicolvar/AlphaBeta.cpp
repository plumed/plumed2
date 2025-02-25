/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2023 The plumed team
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
#include "core/ActionWithValue.h"
#include "core/ActionRegister.h"
#include "core/PlumedMain.h"
#include "core/ActionSet.h"
#include <string>
#include <cmath>

//+PLUMEDOC MCOLVAR ALPHABETA
/*
Calculate the alpha beta CV

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class AlphaBeta : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit AlphaBeta(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(AlphaBeta,"ALPHABETA")

void AlphaBeta::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.add("numbered","ATOMS","the atoms involved for each of the torsions you wish to calculate. "
           "Keywords like ATOMS1, ATOMS2, ATOMS3,... should be listed and one torsion will be "
           "calculated for each ATOM keyword you specify");
  keys.reset_style("ATOMS","atoms");
  keys.add("numbered","REFERENCE","the reference values for each of the torsional angles.  If you use a single REFERENCE value the "
           "same reference value is used for all torsions");
  keys.add("numbered","COEFFICIENT","the coefficient for each of the torsional angles.  If you use a single COEFFICIENT value the "
           "same reference value is used for all torsional angles");
  keys.setValueDescription("scalar","the alpha beta CV");
  keys.needsAction("CONSTANT");
  keys.needsAction("TORSION");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
  keys.needsAction("SUM");
}

AlphaBeta::AlphaBeta(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the reference value
  std::string refstr;
  parse("REFERENCE",refstr);
  unsigned nref=0;
  if( refstr.length()==0 ) {
    for(unsigned i=0;; ++i) {
      std::string refval;
      if( !parseNumbered( "REFERENCE", i+1, refval ) ) {
        break;
      }
      if( i==0 ) {
        refstr = refval;
      } else {
        refstr += "," + refval;
      }
      nref++;
    }
  }
  std::string coeffstr;
  parse("COEFFICIENT",coeffstr);
  unsigned ncoeff=0;
  if( coeffstr.length()==0 ) {
    for(unsigned i=0;; ++i) {
      std::string coeff;
      if( !parseNumbered( "COEFFICIENT", i+1, coeff) ) {
        break;
      }
      if( i==0 ) {
        coeffstr = coeff;
      } else {
        coeffstr += "," + coeff;
      }
      ncoeff++;
    }
  }
  if( coeffstr.length()==0 ) {
    coeffstr="1";
  }
  // Calculate angles
  readInputLine( getShortcutLabel() + "_torsions: TORSION " + convertInputLineToString() );
  ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( getShortcutLabel() + "_torsions" );
  plumed_assert( av && (av->copyOutput(0))->getRank()==1 );
  if( nref==0 ) {
    std::string refval=refstr;
    for(unsigned i=1; i<(av->copyOutput(0))->getShape()[0]; ++i) {
      refstr += "," + refval;
    }
  } else if( nref!=(av->copyOutput(0))->getShape()[0] ) {
    error("mismatch between number of reference values and number of ATOMS specified");
  }
  if( ncoeff==0 ) {
    std::string coeff=coeffstr;
    for(unsigned i=1; i<(av->copyOutput(0))->getShape()[0]; ++i) {
      coeffstr += "," + coeff;
    }
  } else if( ncoeff!=(av->copyOutput(0))->getShape()[0] ) {
    error("mismatch between number of coefficients and number of ATOMS specified");
  }
  readInputLine( getShortcutLabel() + "_ref: CONSTANT VALUES=" + refstr );
  readInputLine( getShortcutLabel() + "_coeff: CONSTANT VALUES=" + coeffstr );
  // Caculate difference from reference using combine
  readInputLine( getShortcutLabel() + "_comb: COMBINE ARG=" + getShortcutLabel() + "_torsions," + getShortcutLabel() + "_ref COEFFICIENTS=1,-1 PERIODIC=NO" );
  // Now matheval for cosine bit
  readInputLine( getShortcutLabel() + "_cos: CUSTOM ARG=" + getShortcutLabel() + "_comb," + getShortcutLabel() + "_coeff FUNC=y*(0.5+0.5*cos(x)) PERIODIC=NO");
  // And combine to get final value
  readInputLine( getShortcutLabel() + ": SUM ARG=" + getShortcutLabel() + "_cos PERIODIC=NO");
}

}
}
