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

//+PLUMEDOC MCOLVAR PAIRENTROPY
/*
Calculate the KL Entropy from the radial distribution function

This shortcut provides an implementation of the CV that is described in the paper that is cited below. Although, the implementation
that is offered here is not the one that was used to produce the results described in that paper the values produced by this
implementation have been tested against the implementation that was used in the paper, which can be found [here](https://sites.google.com/site/pablompiaggi/scripts/pair-entropy/pair-entropy-cv?authuser=0).

An example input that calculates and prints the PAIRENTROPY CV is shown below:

```plumed
pp: PAIRENTROPY ...
   GROUP=1-108 MAXR=2.0 GRID_BIN=20
   BANDWIDTH=0.13 KERNEL=gaussian
   CUTOFF=6.25
...
PRINT ARG=pp FILE=colvar
```

By expanding the shortct in the input above you can see how features that are already avaialble within PLUMED can be reused to calculate this CV.  The resulting implementation is likely slower than the
direct implementation that is available [here](https://sites.google.com/site/pablompiaggi/scripts/pair-entropy/pair-entropy-cv?authuser=0).  We hope, however, that this implementation helps others to understand
how this CV is constructed.

If you would like to calculate the pair entropy based on the radial distribution of the atoms in GROUPB around GROUPA you can use an input similar to the one shown below:

```plumed
pp: PAIRENTROPY ...
   GROUPA=1-108 GROUPB=109-300
   GRID_BIN=20 MAXR=2.0
   BANDWIDTH=0.13 KERNEL=gaussian
   CUTOFF=6.25 DENSITY=1.0
...
PRINT ARG=pp FILE=colvar
```

Notice that we have also used the `DENSITY` keyword to set the background density that is used when normalizing the radial distribution function explicity to 1 atom$/nm^{3}$.
When this keyword is not used, this density is calculated by dividing the number of atoms by the volume of the box as you can see if you expand the shortcut in the
first input above.

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace gridtools {

class PairEntropy : public ActionShortcut {
public:
  static void registerKeywords( Keywords& keys );
  explicit PairEntropy(const ActionOptions&ao);
};

PLUMED_REGISTER_ACTION(PairEntropy,"PAIRENTROPY")

void PairEntropy::registerKeywords( Keywords& keys ) {
  RDF::registerKeywords( keys );
  keys.needsAction("RDF");
  keys.addDeprecatedKeyword("SIGMA","BANDWIDTH");
  keys.setValueDescription("scalar","the KL-entropy that is computed from the radial distribution function");
  keys.addDOI("10.1103/PhysRevLett.119.015701");
  keys.needsAction("INTERPOLATE_GRID");
  keys.needsAction("INTEGRATE_GRID");
}

PairEntropy::PairEntropy(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string ref_str, ref_name;
  parse("REFERENCE",ref_name);
  if( ref_name.length()>0 ) {
    ref_str = "REFERENCE=" + ref_name;
  } else {
    ref_name = getShortcutLabel() + "_rdf";
  }
  // Read in the atoms and get the number of atoms that we are using
  std::string atom_str, group_str, natoms;
  parse("GROUP",group_str);
  if( group_str.length()>0 ) {
    atom_str="GROUP=" + group_str;
    std::vector<std::string> awords=Tools::getWords(group_str,"\t\n ,");
    Tools::interpretRanges( awords );
    Tools::convert( awords.size(), natoms );
  } else {
    std::string groupa_str, groupb_str;
    parse("GROUPA",groupa_str);
    parse("GROUPB",groupb_str);
    atom_str="GROUPA=" + groupa_str + " GROUPB=" + groupb_str;
    std::vector<std::string> awords=Tools::getWords(groupb_str,"\t\n ,");
    Tools::interpretRanges( awords );
    Tools::convert( awords.size()+1, natoms );
  }
  // Read in all other keywords and create the RDF object
  std::string maxr, nbins, dens, bw="", cutoff, kernel;
  parse("MAXR",maxr);
  parse("GRID_BIN",nbins);
  parse("DENSITY",dens);
  parse("BANDWIDTH",bw);
  if( bw.length()==0 ) {
    parse("SIGMA",bw);
  }
  parse("CUTOFF",cutoff);
  parse("KERNEL",kernel);
  std::string dens_str;
  if( dens.length()>0 ) {
    dens_str = " DENSITY=" + dens;
  }
  readInputLine( getShortcutLabel() + "_rdf: RDF " + atom_str + " KERNEL=" + kernel + " CUTOFF=" + cutoff + " GRID_BIN=" + nbins + " MAXR=" + maxr + dens_str + " BANDWIDTH=" + bw + " " + ref_str);
  // And compute the two functions we are integrating (we use two matheval objects here and sum them in order to avoid nans from taking logarithms of zero)
  readInputLine( getShortcutLabel() + "_conv_t1: CUSTOM ARG=" + getShortcutLabel() + "_rdf," + ref_name + "_x2 FUNC=x*y*log(x) PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_conv_t2: CUSTOM ARG=" + getShortcutLabel() + "_rdf," + ref_name + "_x2 FUNC=(1-x)*y PERIODIC=NO");
  readInputLine( getShortcutLabel() + "_conv: CUSTOM ARG=" + getShortcutLabel() + "_conv_t1," + getShortcutLabel() + "_conv_t2 FUNC=x+y PERIODIC=NO");
  // Now integrate using trapezium rule
  readInputLine( getShortcutLabel() + "_midp: INTERPOLATE_GRID ARG=" + getShortcutLabel() + "_conv INTERPOLATION_TYPE=linear MIDPOINTS"); // First interpolate onto midpoints
  readInputLine( getShortcutLabel() + "_int: INTEGRATE_GRID ARG=" + getShortcutLabel() + "_midp PERIODIC=NO"); // And then integrate
  // And multiply by final normalizing constant
  std::string norm_str;
  if( dens.length()>0 ) {
    norm_str = " FUNC=-2*pi*x*" + dens;
  } else {
    norm_str = "," + ref_name + "_vol FUNC=-(2*pi*x/y)*" + natoms;
  }
  readInputLine( getShortcutLabel() + ": CUSTOM PERIODIC=NO ARG=" + getShortcutLabel() + "_int" + norm_str );
}

}
}
