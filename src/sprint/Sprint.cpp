/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2013-2020 The plumed team
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

//+PLUMEDOC MATRIXF SPRINT
/*
Calculate SPRINT topological variables from an adjacency matrix.

The SPRINT topological variables are calculated from the largest eigenvalue, \f$\lambda\f$ of
an \f$n\times n\f$ adjacency matrix and its corresponding eigenvector, \f$\mathbf{V}\f$, using:

\f[
s_i = \sqrt{n} \lambda v_i
\f]

You can use different quantities to measure whether or not two given atoms/molecules are
adjacent or not in the adjacency matrix.  The simplest measure of adjacency is is whether
two atoms/molecules are within some cutoff of each other.  Further complexity can be added by
insisting that two molecules are adjacent if they are within a certain distance of each
other and if they have similar orientations.

\par Examples

This example input calculates the 7 SPRINT coordinates for a 7 atom cluster of Lennard-Jones
atoms and prints their values to a file.  In this input the SPRINT coordinates are calculated
in the manner described in ?? so two atoms are adjacent if they are within a cutoff:

\plumedfile
DENSITY SPECIES=1-7 LABEL=d1
CONTACT_MATRIX ATOMS=d1 SWITCH={RATIONAL R_0=0.1} LABEL=mat
SPRINT MATRIX=mat LABEL=ss
PRINT ARG=ss.* FILE=colvar
\endplumedfile

This example input calculates the 14 SPRINT coordinates for a molecule composed of 7 hydrogen and
7 carbon atoms.  Once again two atoms are adjacent if they are within a cutoff:

\plumedfile
DENSITY SPECIES=1-7 LABEL=c
DENSITY SPECIES=8-14 LABEL=h

CONTACT_MATRIX ...
  ATOMS=c,h
  SWITCH11={RATIONAL R_0=2.6 NN=6 MM=12}
  SWITCH12={RATIONAL R_0=2.2 NN=6 MM=12}
  SWITCH22={RATIONAL R_0=2.2 NN=6 MM=12}
  LABEL=mat
... CONTACT_MATRIX

SPRINT MATRIX=mat LABEL=ss

PRINT ARG=ss.* FILE=colvar
\endplumedfile

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace sprint {

class Sprint : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Sprint(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Sprint,"SPRINT")

void Sprint::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("optional","MATRIX","the matrix that you would like to perform SPRINT on");
  keys.add("numbered","GROUP","specifies the list of atoms that should be assumed indistinguishable");
  keys.add("numbered","SWITCH","specify the switching function to use between two sets of indistinguishable atoms");
  keys.needsAction("CONTACT_MATRIX"); keys.needsAction("DIAGONALIZE"); keys.needsAction("CUSTOM");
  keys.needsAction("SELECT_COMPONENTS"); keys.needsAction("SORT"); keys.needsAction("COMBINE");
  keys.addOutputComponent("coord","default","the sprint coordinates");
}

Sprint::Sprint(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao)
{
  std::string matinp; parse("MATRIX",matinp);
  if( matinp.length()==0 ) {
    readInputLine( getShortcutLabel() + "_jmat: CONTACT_MATRIX " + convertInputLineToString() );
    matinp = getShortcutLabel() + "_jmat";
  }
  std::vector<unsigned> nin_group; unsigned ntot_atoms=0;
  for(unsigned i=1;; ++i) {
    std::string inum; Tools::convert( i, inum );
    ActionWithValue* av = plumed.getActionSet().selectWithLabel<ActionWithValue*>( matinp + inum + inum );
    if( !av ) break ;
    unsigned natoms = (av->copyOutput(0))->getShape()[0]; nin_group.push_back( natoms ); ntot_atoms += natoms;
  }

  // Diagonalization
  readInputLine( getShortcutLabel() + "_diag: DIAGONALIZE ARG=" + matinp + " VECTORS=1");
  // Compute sprint coordinates as product of eigenvalue and eigenvector times square root of number of atoms in all groups
  std::string str_natoms; Tools::convert( ntot_atoms, str_natoms );
  readInputLine( getShortcutLabel() + "_sp: CUSTOM ARG=" + getShortcutLabel() + "_diag.vals-1," + getShortcutLabel() +
                 "_diag.vecs-1 FUNC=sqrt(" + str_natoms + ")*x*y PERIODIC=NO");
  // Sort sprint coordinates for each group of atoms
  unsigned k=0, kk=0;
  for(unsigned j=0; j<nin_group.size(); ++j) {
    std::string jnum, knum; Tools::convert( j+1, jnum ); Tools::convert(k+1, knum); k++;
    std::string sort_act = getShortcutLabel() + "_selection" + jnum + ": SELECT_COMPONENTS ARG=" + getShortcutLabel() + "_sp COMPONENTS=" + knum;
    for(unsigned n=1; n<nin_group[j]; ++n) {
      Tools::convert( k+1, knum ); sort_act += ","+ knum; k++;
    }
    readInputLine( sort_act );
    readInputLine( getShortcutLabel() + jnum + ": SORT ARG=" + getShortcutLabel() + "_selection" + jnum );
    for(unsigned n=0; n<nin_group[j]; ++n) {
      std::string knum, nnum; Tools::convert( kk, knum ); Tools::convert( n+1, nnum ); kk++;
      readInputLine( getShortcutLabel() + "_coord-" + knum + ": COMBINE ARG=" + getShortcutLabel() + jnum + "." + nnum + " PERIODIC=NO" );
    }
  }
}

}
}
