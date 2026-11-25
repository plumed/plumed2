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
Calculate an angle.

This command can be used to compute the angle between three atoms. Alternatively
if four atoms appear in the atom
specification it calculates the angle between
two vectors identified by two pairs of atoms.

If _three_ atoms are given, the angle is defined as:
\f[
\theta=\arccos\left(\frac{ {\bf r}_{21}\cdot {\bf r}_{23}}{
|{\bf r}_{21}| |{\bf r}_{23}|}\right)
\f]
Here \f$ {\bf r}_{ij}\f$ is the distance vector among the
i-th and the j-th listed atom.

If _four_ atoms are given, the angle is defined as:
\f[
\theta=\arccos\left(\frac{ {\bf r}_{21}\cdot {\bf r}_{34}}{
|{\bf r}_{21}| |{\bf r}_{34}|}\right)
\f]

Notice that angles defined in this way are non-periodic variables and
their value is limited by definition between 0 and \f$\pi\f$.

The vectors \f$ {\bf r}_{ij}\f$ are by default evaluated taking
periodic boundary conditions into account.
This behavior can be changed with the NOPBC flag.

\par Examples

This command tells plumed to calculate the angle between the vector connecting atom 1 to atom 2 and
the vector connecting atom 2 to atom 3 and to print it on file COLVAR1. At the same time,
the angle between vector connecting atom 1 to atom 2 and the vector connecting atom 3 to atom 4 is printed
on file COLVAR2.
\plumedfile

a: ANGLE ATOMS=1,2,3
# equivalently one could state:
# a: ANGLE ATOMS=1,2,2,3

b: ANGLE ATOMS=1,2,3,4

PRINT ARG=a FILE=COLVAR1
PRINT ARG=b FILE=COLVAR2
\endplumedfile

This final example instructs plumed to calculate all the angles in the first coordination
spheres of the atoms. The bins for a normalized histogram of the distribution is then output

\plumedfile
ANGLES GROUP=1-38 HISTOGRAM={GAUSSIAN LOWER=0.0 UPPER=pi NBINS=20} SWITCH={GAUSSIAN R_0=1.0} LABEL=a1
PRINT ARG=a1.* FILE=colvar
\endplumedfile


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
  MultiColvarShortcuts::shortcutKeywords( keys );
  keys.needsAction("ANGLE");
  keys.needsAction("COORD_ANGLES");
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



