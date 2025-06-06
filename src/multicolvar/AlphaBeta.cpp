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
Measures a distance including pbc between the instantaneous values of a set of torsional angles and set of reference values.

This shortcut calculates the following quantity.

$$
s = \frac{1}{2} \sum_i c_i \left[ 1 + \cos( \phi_i - \phi_i^{\textrm{Ref}} ) \right]
$$

where the $\phi_i$ values are the instantaneous values for the [TORSION](TORSION.md) angles of interest.
The $\phi_i^{\textrm{Ref}}$ values are reference values for the torsional angles that are specified in the input file.

The following provides an example of the input for an alpha beta similarity.

```plumed
ab: ALPHABETA ...
ATOMS1=168,170,172,188 REFERENCE1=3.14
ATOMS2=170,172,188,190 REFERENCE2=3.14
ATOMS3=188,190,192,230 REFERENCE3=3.14
...
PRINT ARG=ab FILE=colvar STRIDE=10
```

Because all the reference values are the same we can also calculate the same quantity using

```plumed
ab: ALPHABETA ...
ATOMS1=168,170,172,188 REFERENCE=3.14
ATOMS2=170,172,188,190
ATOMS3=188,190,192,230
...
PRINT ARG=ab FILE=colvar STRIDE=10
```

Writing out the atoms involved in all the torsion angles in this way can be rather tedious. Thankfully if you are working with protein you
can avoid this by using the [MOLINFO](MOLINFO.md) command.  PLUMED uses the pdb file that you provide to this command to learn
about the topology of the protein molecule.  This means that you can specify torsion angles using the following syntax:

```plumed
#SETTINGS MOLFILE=regtest/basic/rt32/helix.pdb
MOLINFO MOLTYPE=protein STRUCTURE=regtest/basic/rt32/helix.pdb
ab: ALPHABETA ...
ATOMS1=@phi-3 REFERENCE=3.14 COEFFICIENT1=2
ATOMS2=@psi-3                COEFFICIENT2=0.5
ATOMS3=@phi-4                COEFFICIENT3=1
...
PRINT ARG=ab FILE=colvar STRIDE=10
```

Here, `@phi-3` tells plumed that you would like to calculate the $\phi$ angle in the third residue of the protein.
Similarly `@psi-4` tells plumed that you want to calculate the $\psi$ angle of the fourth residue of the protein.
Notice, also, that in the first two examples the coefficients $c_i$ in the expression above were all set equal to one.
In the example above we use the COEFFICIENT keywords to set these quantities to three different values.

Notice, last of all, that in the above examples we reassemble any molecules that have been broken by the periodic boundary
conditions using a procedure like that used in [WHOLEMOLECULES](WHOLEMOLECULES.md) before calculating the torsion angles.
If you wish to turn this off for any reason you use the NOPBC flag as shown below:

```plumed
ab: ALPHABETA ...
ATOMS1=168,170,172,188 REFERENCE=3.14
ATOMS2=170,172,188,190 NOPBC
ATOMS3=188,190,192,230
...
PRINT ARG=ab FILE=colvar STRIDE=10
```


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
