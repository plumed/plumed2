/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2023 The plumed team
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

//+PLUMEDOC MCOLVAR DISTANCES
/*
Calculate the distances between multiple piars of atoms

\par Examples

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR XDISTANCES
/*
Calculate the x components of the vectors connecting one or many pairs of atoms.
You can then calculate functions of the distribution of values such as the minimum, the number less than a certain quantity and so on.

\par Examples

The following input tells plumed to calculate the x-component of the vector connecting atom 3 to atom 5 and
the x-component of the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then
printed
\plumedfile
d1: XDISTANCES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min
\endplumedfile
(See also \ref PRINT).


The following input tells plumed to calculate the x-component of the vector connecting atom 3 to atom 5 and
the x-component of the vector connecting atom 1 to atom 2.  The number of values that are
less than 0.1nm is then printed to a file.
\plumedfile
d1: XDISTANCES ATOMS1=3,5 ATOMS2=1,2 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan
\endplumedfile
(See also \ref PRINT \ref switchingfunction).

The following input tells plumed to calculate the x-components of all the distinct vectors that can be created
between atoms 1, 2 and 3 (i.e. the vectors between atoms 1 and 2, atoms 1 and 3 and atoms 2 and 3).
The average of these quantities is then calculated.
\plumedfile
d1: XDISTANCES GROUP=1-3 MEAN
PRINT ARG=d1.mean
\endplumedfile
(See also \ref PRINT)

The following input tells plumed to calculate all the vectors connecting the the atoms in GROUPA to the atoms in GROUPB.
In other words the vector between atoms 1 and 2 and the vector between atoms 1 and 3.  The number of values
more than 0.1 is then printed to a file.
\plumedfile
d1: XDISTANCES GROUPA=1 GROUPB=2,3 MORE_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.morethan
\endplumedfile
(See also \ref PRINT \ref switchingfunction)
*/
//+ENDPLUMEDOC


//+PLUMEDOC MCOLVAR YDISTANCES
/*
Calculate the y components of the vectors connecting one or many pairs of atoms.
You can then calculate functions of the distribution of values such as the minimum, the number less than a certain quantity and so on.

\par Examples

The following input tells plumed to calculate the y-component of the vector connecting atom 3 to atom 5 and
the y-component of the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then
printed
\plumedfile
d1: YDISTANCES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min
\endplumedfile
(See also \ref PRINT).


The following input tells plumed to calculate the y-component of the vector connecting atom 3 to atom 5 and
the y-component of the vector connecting atom 1 to atom 2.  The number of values that are
less than 0.1nm is then printed to a file.
\plumedfile
d1: YDISTANCES ATOMS1=3,5 ATOMS2=1,2 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan
\endplumedfile
(See also \ref PRINT \ref switchingfunction).

The following input tells plumed to calculate the y-components of all the distinct vectors that can be created
between atoms 1, 2 and 3 (i.e. the vectors between atoms 1 and 2, atoms 1 and 3 and atoms 2 and 3).
The average of these quantities is then calculated.
\plumedfile
d1: YDISTANCES GROUP=1-3 MEAN
PRINT ARG=d1.mean
\endplumedfile
(See also \ref PRINT)

The following input tells plumed to calculate all the vectors connecting the the atoms in GROUPA to the atoms in GROUPB.
In other words the vector between atoms 1 and 2 and the vector between atoms 1 and 3.  The number of values
more than 0.1 is then printed to a file.
\plumedfile
d1: YDISTANCES GROUPA=1 GROUPB=2,3 MORE_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.morethan
\endplumedfile
(See also \ref PRINT \ref switchingfunction)

*/
//+ENDPLUMEDOC

//+PLUMEDOC MCOLVAR ZDISTANCES
/*
Calculate the z components of the vectors connecting one or many pairs of atoms.
You can then calculate functions of the distribution of values such as the minimum, the number less than a certain quantity and so on.

\par Examples

The following input tells plumed to calculate the z-component of the vector connecting atom 3 to atom 5 and
the z-component of the vector connecting atom 1 to atom 2.  The minimum of these two quantities is then
printed
\plumedfile
d1: ZDISTANCES ATOMS1=3,5 ATOMS2=1,2 MIN={BETA=0.1}
PRINT ARG=d1.min
\endplumedfile
(See also \ref PRINT).


The following input tells plumed to calculate the z-component of the vector connecting atom 3 to atom 5 and
the z-component of the vector connecting atom 1 to atom 2.  The number of values that are
less than 0.1nm is then printed to a file.
\plumedfile
d1: ZDISTANCES ATOMS1=3,5 ATOMS2=1,2 LESS_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.lessthan
\endplumedfile
(See also \ref PRINT \ref switchingfunction).

The following input tells plumed to calculate the z-components of all the distinct vectors that can be created
between atoms 1, 2 and 3 (i.e. the vectors between atoms 1 and 2, atoms 1 and 3 and atoms 2 and 3).
The average of these quantities is then calculated.
\plumedfile
d1: ZDISTANCES GROUP=1-3 MEAN
PRINT ARG=d1.mean
\endplumedfile
(See also \ref PRINT)

The following input tells plumed to calculate all the vectors connecting the the atoms in GROUPA to the atoms in GROUPB.
In other words the vector between atoms 1 and 2 and the vector between atoms 1 and 3.  The number of values
more than 0.1 is then printed to a file.
\plumedfile
d1: ZDISTANCES GROUPA=1 GROUPB=2,3 MORE_THAN={RATIONAL R_0=0.1}
PRINT ARG=d1.morethan
\endplumedfile
(See also \ref PRINT \ref switchingfunction)

*/
//+ENDPLUMEDOC


namespace PLMD {
namespace multicolvar {

class Distances : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Distances(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Distances,"DISTANCES")
PLUMED_REGISTER_ACTION(Distances,"XDISTANCES")
PLUMED_REGISTER_ACTION(Distances,"YDISTANCES")
PLUMED_REGISTER_ACTION(Distances,"ZDISTANCES")

void Distances::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("atoms-1","GROUP","Calculate the distance between each distinct pair of atoms in the group");
  keys.add("atoms-2","GROUPA","Calculate the distances between all the atoms in GROUPA and all "
           "the atoms in GROUPB. This must be used in conjunction with GROUPB.");
  keys.add("atoms-2","GROUPB","Calculate the distances between all the atoms in GROUPA and all the atoms "
           "in GROUPB. This must be used in conjunction with GROUPA.");
  keys.add("numbered","ATOMS","the pairs of atoms that you would like to calculate the angles for");
  keys.addFlag("NOPBC",false,"ignore the periodic boundary conditions when calculating distances");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the distance separately and store them as label.x, label.y and label.z");
  keys.addFlag("SCALED_COMPONENTS",false,"calculate the a, b and c scaled components of the distance separately and store them as label.a, label.b and label.c");
  keys.addFlag("LOWMEM",false,"this flag does nothing and is present only to ensure back-compatibility");
  keys.reset_style("ATOMS","atoms");
  MultiColvarShortcuts::shortcutKeywords( keys );
  keys.add("atoms","ORIGIN","calculate the distance of all the atoms specified using the ATOMS keyword from this point");
  keys.add("numbered","LOCATION","the location at which the CV is assumed to be in space");
  keys.reset_style("LOCATION","atoms");
  keys.setValueDescription("vector","the DISTANCES between the each pair of atoms that were specified");
  keys.addOutputComponent("x","COMPONENTS","vector","the x-components of the distance vectors");
  keys.addOutputComponent("y","COMPONENTS","vector","the y-components of the distance vectors");
  keys.addOutputComponent("z","COMPONENTS","vector","the z-components of the distance vectors");
  keys.needsAction("GROUP");
  keys.needsAction("DISTANCE");
  keys.needsAction("CENTER");
}

Distances::Distances(const ActionOptions& ao):
  Action(ao),
  ActionShortcut(ao) {
  // Create distances
  bool lowmem;
  parseFlag("LOWMEM",lowmem);
  if( lowmem ) {
    warning("LOWMEM flag is deprecated and is no longer required for this action");
  }
  std::string dline = getShortcutLabel() + ": DISTANCE";
  bool nopbc;
  parseFlag("NOPBC",nopbc);
  if( nopbc ) {
    dline += " NOPBC";
  }
  if( getName()=="DISTANCES" ) {
    bool comp;
    parseFlag("COMPONENTS",comp);
    if( comp ) {
      dline += " COMPONENTS";
    }
    bool scomp;
    parseFlag("SCALED_COMPONENTS",scomp);
    if( scomp ) {
      dline += " SCALED_COMPONENTS";
    }
  } else {
    dline += " COMPONENTS";
  }
  // Parse origin
  std::string num, ostr;
  parse("ORIGIN",ostr);
  if( ostr.length()>0 ) {
    // Parse atoms
    std::vector<std::string> afstr;
    MultiColvarShortcuts::parseAtomList("ATOMS",afstr,this);
    for(unsigned i=0; i<afstr.size(); ++i) {
      Tools::convert( i+1, num );
      dline += " ATOMS" + num + "=" + ostr + "," + afstr[i];
    }
  } else {
    std::vector<std::string> grp;
    MultiColvarShortcuts::parseAtomList("GROUP",grp,this);
    std::vector<std::string> grpa;
    MultiColvarShortcuts::parseAtomList("GROUPA",grpa,this);
    if( grp.size()>0 ) {
      if( grpa.size()>0 ) {
        error("should not be using GROUPA in tandem with GROUP");
      }
      unsigned n=0;
      for(unsigned i=1; i<grp.size(); ++i) {
        for(unsigned j=0; j<i; ++j) {
          std::string num;
          Tools::convert( n+1, num );
          n++;
          dline += " ATOMS" + num + "=" + grp[i] + "," + grp[j];
        }
      }
    } else if( grpa.size()>0 ) {
      std::vector<std::string> grpb;
      MultiColvarShortcuts::parseAtomList("GROUPB",grpb,this);
      if( grpb.size()==0 ) {
        error("found GROUPA but no corresponding GROUPB");
      }
      std::string grpstr = getShortcutLabel() + "_grp: GROUP ATOMS=";
      bool printcomment=false;
      for(unsigned i=0; i<grpa.size(); ++i) {
        for(unsigned j=0; j<grpb.size(); ++j) {
          std::string num;
          Tools::convert( i*grpb.size() + j + 1, num );
          dline += " ATOMS" + num + "=" + grpa[i] + "," + grpb[j];
          if( i*grpb.size() + j<6 ) {
            readInputLine( getShortcutLabel() + "_vatom" + num + ": CENTER ATOMS=" + grpa[i] + "," + grpb[j], true );
          } else {
            readInputLine( getShortcutLabel() + "_vatom" + num + ": CENTER ATOMS=" + grpa[i] + "," + grpb[j], false );
            printcomment=true;
          }
          if( i+j==0 ) {
            grpstr += getShortcutLabel() + "_vatom" + num;
          } else {
            grpstr += "," + getShortcutLabel() + "_vatom" + num;
          }
        }
      }
      std::string num;
      Tools::convert( grpa.size()*grpb.size(), num );
      if( printcomment ) {
        addCommentToShortcutOutput("# A further " + num + " CENTER like the ones above were also created but are not shown");
      }
      readInputLine( grpstr );
    } else {
      std::string grpstr = getShortcutLabel() + "_grp: GROUP ATOMS=";
      for(unsigned i=1;; ++i) {
        std::string atstring;
        parseNumbered("ATOMS",i,atstring);
        if( atstring.length()==0 ) {
          break;
        }
        std::string locstr;
        parseNumbered("LOCATION",i,locstr);
        if( locstr.length()==0 ) {
          std::string num;
          Tools::convert( i, num );
          readInputLine( getShortcutLabel() + "_vatom" + num + ": CENTER ATOMS=" + atstring );
          if( i==1 ) {
            grpstr += getShortcutLabel() + "_vatom" + num;
          } else {
            grpstr += "," + getShortcutLabel() + "_vatom" + num;
          }
        } else {
          if( i==1 ) {
            grpstr += locstr;
          } else {
            grpstr += "," + locstr;
          }
        }
        std::string num;
        Tools::convert( i, num );
        dline += " ATOMS" + num + "=" + atstring;
      }
      readInputLine( grpstr );
    }
  }
  readInputLine( dline );
  // Add shortcuts to label
  if( getName()=="DISTANCES" ) {
    MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
  } else if( getName()=="XDISTANCES" ) {
    MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + ".x", "", this );
  } else if( getName()=="YDISTANCES" ) {
    MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + ".y", "", this );
  } else if( getName()=="ZDISTANCES" ) {
    MultiColvarShortcuts::expandFunctions( getShortcutLabel(), getShortcutLabel() + ".z", "", this );
  }
}

}
}
