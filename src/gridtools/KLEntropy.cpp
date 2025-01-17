/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2016-2023 The plumed team
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
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"

//+PLUMEDOC ANALYSIS KL_ENTROPY
/*
Calculate the KL entropy of a distribution

\par Examples


*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class KLEntropy : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit KLEntropy(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(KLEntropy,"KL_ENTROPY")

void KLEntropy::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords(keys);
  keys.add("compulsory","ARG","the grid that you wish to use in the KL entropy calculation");
  keys.add("compulsory","REFERENCE","a file containing the reference density in grid format");
  keys.add("compulsory","VALUE","the name of the value that should be read from the grid");
  keys.setValueDescription("the value of the KL-Entropy");
  keys.needsAction("REFERENCE_GRID");
  keys.needsAction("CUSTOM");
  keys.needsAction("INTEGRATE_GRID");
}

KLEntropy::KLEntropy( const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Reference grid object
  std::string ref_str, val_str, input_g;
  parse("VALUE",val_str);
  parse("REFERENCE",ref_str);
  parse("ARG",input_g);
  readInputLine( getShortcutLabel() + "_ref: REFERENCE_GRID VALUE=" + val_str + " FILE=" + ref_str );
  // Compute KL divergence
  if( input_g=="") {
    plumed_merror("could not find ARG keyword in input to KL_ENTROPY");
  }
  readInputLine( getShortcutLabel()  + "_kl: CUSTOM ARG=" + input_g + "," + getShortcutLabel() + "_ref FUNC=y*log(y/(0.5*(x+y))) PERIODIC=NO");
  // Compute integral
  readInputLine( getShortcutLabel() + ": INTEGRATE_GRID ARG=" + getShortcutLabel() + "_kl PERIODIC=NO");
}

}
}
