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
#include "RDF.h"
#include "core/ActionRegister.h"
#include "core/ActionShortcut.h"

//+PLUMEDOC MCOLVAR PAIRENTROPIES
/*
Calculate the KL entropy from the RDF around each of the atoms

This shortcut can be used to calculate the descriptors of local structure that are described in the paper that is cited below.  This is not the implementation
that was used when that paper was written.  However, the values that this implementation generates were tested against the values that the implementation that
was used in the paper, which can be found [here](https://sites.google.com/site/pablompiaggi/scripts/pair-entropy/pair-entropy-fingerprint?authuser=0).

The following example illustrates how this shortcut can be used to calculate and print the descriptor for a set of 64 atoms:

```plumed
pp: PAIRENTROPIES ATOMS=1-64 MAXR=2 GRID_BIN=20 BANDWIDTH=0.13
DUMPATOMS ATOMS=1-64 ARG=pp FILE=p_entropies.xyz
```

If you expand the shortcut you will notice that this shortcut creates __a lot__ of actions.  It is thus likely to be considerably slower than the implementation
[here](https://sites.google.com/site/pablompiaggi/scripts/pair-entropy/pair-entropy-fingerprint?authuser=0). However, we hope that the implementation here is a
useful reference implementation that anyone interested in using this method can use to test their implementations of it. In addition, we hope that the shortcut
illustrates how the calculation of this complex descriptor can be broken down into a set of simpler calculations.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class PairEntropies : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit PairEntropies(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(PairEntropies,"PAIRENTROPIES")

void PairEntropies::registerKeywords( Keywords& keys ) {
  RDF::registerKeywords( keys );
  keys.remove("GROUP");
  keys.remove("GROUPA");
  keys.remove("GROUPB");
  keys.remove("ATOMS");
  keys.add("atoms","ATOMS","the atoms that you would like to compute the entropies for");
  keys.setValueDescription("vector","the a vector containing the KL-entropy that is computed from the radial distribution function around each of the atoms in the input");
  keys.addDOI("10.1063/1.4998408");
  keys.needsAction("PAIRENTROPY");
  keys.needsAction("INTERPOLATE_GRID");
  keys.needsAction("INTEGRATE_GRID");
  keys.needsAction("CUSTOM");
  keys.needsAction("CONCATENATE");
}

PairEntropies::PairEntropies(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  // Read in the atoms involved
  std::string atoms_str;
  parse("ATOMS",atoms_str);
  // Create the x2 value
  std::string maxr, nbins, dens;
  parse("MAXR",maxr);
  parse("GRID_BIN",nbins);
  parse("DENSITY",dens);
  std::string grid_setup = "GRID_MIN=0 GRID_MAX=" + maxr + " GRID_BIN=" + nbins;
  RDF::createX2ReferenceObject( getShortcutLabel(), grid_setup, dens.length()==0, this );
  // Create the number of input atoms
  std::string pair_str = "GRID_BIN=" + nbins + " MAXR=" + maxr + " " + convertInputLineToString();
  std::vector<std::string> awords=Tools::getWords(atoms_str,"\t\n ,");
  Tools::interpretRanges( awords );
  // Now create all the pairentropy objects
  for(unsigned i=0; i<awords.size(); ++i) {
    std::string ilab, jlab, atoms_str2;
    Tools::convert( awords[i], ilab );
    for(unsigned j=0; j<awords.size(); ++j) {
      if( awords[i]==awords[j] ) {
        continue;
      }
      Tools::convert( awords[j], jlab );
      if( atoms_str2.length()==0 ) {
        atoms_str2 = jlab;
      } else {
        atoms_str2 += "," + jlab;
      }
    }
    readInputLine( getShortcutLabel() + ilab + ": PAIRENTROPY GROUPA=" + ilab + " GROUPB=" + atoms_str2 + " " + pair_str + " REFERENCE=" + getShortcutLabel() );
  }
  // And compose a vector containing the values of all the pair entropies
  std::string num, argstr = "ARG=" + getShortcutLabel() + "1";
  for(unsigned i=1; i<awords.size(); ++i) {
    Tools::convert( i+1, num );
    argstr += "," + getShortcutLabel() + num;
  }
  readInputLine( getShortcutLabel() + ": CONCATENATE " + argstr );
}

}
}
