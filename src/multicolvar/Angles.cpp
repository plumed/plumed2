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
#include "core/ActionRegister.h"
#include "MultiColvarShortcuts.h"

#include <string>
#include <cmath>

namespace PLMD {
namespace multicolvar {

//+PLUMEDOC COLVAR ANGLES
/*
Calculate functions of the distribution of angles.

__This shortcut action allows you to calculate function of a distribution of angles and reproduces the syntax in older PLUMED versions.
If you look at the example inputs below you can
see how the new syntax operates. We would strongly encourage you to use the newer syntax as it offers greater flexibility.__

You can use this command to calculate functions, $g$, such as:

$$
 f(x) = \sum_{ijk} g( \theta_{ijk} )
$$

where $\theta_{ijk}$ is the angle between the vector connecting atom $i$ and and $j$ and the vector connecting atom $j$ and atom $k$.
Alternatively you can use this command to calculate functions such as:

$$
f(x) = \sum_{ijk} s(r_{ij})s(r_{jk}) g(\theta_{ijk})
$$

where $s(r)$ is a switching function.  This second form means that you can
use this to calculate functions of the angles in the first coordination sphere of
an atom / molecule and hence use the CVs described in the paper in the bibliography below.

The following example tells plumed to calculate all angles involving
at least one atom from GROUPA and two atoms from GROUPB in which the distances
are less than 1.0. The number of angles between $\frac{\pi}{4}$ and
$\frac{3\pi}{4}$ is then output

```plumed
a1: ANGLES GROUPA=1-10 GROUPB=11-100 BETWEEN={GAUSSIAN LOWER=0.25pi UPPER=0.75pi} SWITCH={GAUSSIAN R_0=1.0}
PRINT ARG=a1.between FILE=colvar
```

*/
//+ENDPLUMEDOC

class Angles : public ActionShortcut {
public:
  explicit Angles(const ActionOptions&);
  static void registerKeywords( Keywords& keys );
};

PLUMED_REGISTER_ACTION(Angles,"ANGLES")

void Angles::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms-1","GROUP","Calculate angles for each distinct set of three atoms in the group");
  keys.add("atoms-2","GROUPA","A group of central atoms about which angles should be calculated");
  keys.add("atoms-2","GROUPB","When used in conjunction with GROUPA this keyword instructs plumed "
           "to calculate all distinct angles involving one atom from GROUPA "
           "and two atoms from GROUPB. The atom from GROUPA is the central atom.");
  keys.add("atoms-3","GROUPC","This must be used in conjunction with GROUPA and GROUPB.  All angles "
           "involving one atom from GROUPA, one atom from GROUPB and one atom from "
           "GROUPC are calculated. The GROUPA atoms are assumed to be the central "
           "atoms");
  keys.add("optional","SWITCH","the switching function specifies that only those bonds that have a length that is less than a certain threshold are considered");
  keys.addDOI("https://doi.org/10.1063/1.3628676");
  keys.setValueDescription("vector","the ANGLE for each set of three atoms that were specified");
  MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("ANGLE");
  keys.needsAction("COORD_ANGLES");
  keys.setDeprecated("ANGLE");
}

Angles::Angles(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  std::string swit;
  parse("SWITCH",swit);
  if( swit.length()>0 ) {
    std::string cat, grp;
    parse("GROUPA",cat);
    parse("GROUPB",grp);
    if( cat.length()==0 || grp.length()==0 ) {
      error("must use GROUPA/GROUPB when using SWITCH");
    }
    readInputLine( getShortcutLabel() + ": COORD_ANGLES SWITCH={" +  swit + "} CATOMS=" + cat + " GROUP=" + grp + " " + convertInputLineToString() );
    return;
  }
  std::vector<std::string> group;
  parseVector("GROUP",group);
  std::vector<std::string> groupa;
  parseVector("GROUPA",groupa);
  std::vector<std::string> groupb;
  parseVector("GROUPB",groupb);
  std::vector<std::string> groupc;
  parseVector("GROUPC",groupc);
  if( group.size()>0 ) {
    if( groupa.size()>0 || groupb.size()>0 || groupc.size()>0 ) {
      error("should only be GROUP keyword in input not GROUPA/GROUPB/GROUPC");
    }
    Tools::interpretRanges( group );
    std::string ainput = getShortcutLabel() + ": ANGLE";
    unsigned n=1;
    // Not sure if this triple sum makes any sense
    for(unsigned i=2; i<group.size(); ++i ) {
      for(unsigned j=1; j<i; ++j ) {
        for(unsigned k=0; k<j; ++k) {
          std::string str_n;
          Tools::convert( n, str_n );
          ainput += " ATOMS" + str_n + "=" + group[i] + "," + group[j] + "," + group[k];
          n++;
        }
      }
    }
    readInputLine( ainput );
  } else if( groupc.size()>0 ) {
    Tools::interpretRanges( groupa );
    Tools::interpretRanges( groupb );
    Tools::interpretRanges( groupc );
    unsigned n=1;
    std::string ainput = getShortcutLabel() + ": ANGLE";
    for(unsigned i=0; i<groupa.size(); ++i ) {
      for(unsigned j=0; j<groupb.size(); ++j ) {
        for(unsigned k=0; k<groupc.size(); ++k) {
          std::string str_n;
          Tools::convert( n, str_n );
          ainput += " ATOMS" + str_n + "=" + groupb[j] + "," + groupa[i] + "," + groupc[k];
          n++;
        }
      }
    }
    readInputLine( ainput );
  } else if( groupa.size()>0 ) {
    Tools::interpretRanges( groupa );
    Tools::interpretRanges( groupb );
    unsigned n=1;
    std::string ainput;
    ainput = getShortcutLabel() + ": ANGLE";
    for(unsigned i=0; i<groupa.size(); ++i ) {
      for(unsigned j=1; j<groupb.size(); ++j ) {
        for(unsigned k=0; k<j; ++k) {
          std::string str_n;
          Tools::convert( n, str_n );
          ainput += " ATOMS" + str_n + "=" + groupb[j] + "," + groupa[i] + "," + groupb[k];
          n++;
        }
      }
    }
    readInputLine( ainput );
  }
  MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
}

}
}



