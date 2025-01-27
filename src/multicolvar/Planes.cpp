/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
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
#include "MultiColvarShortcuts.h"
#include "core/ActionShortcut.h"
#include "core/ActionRegister.h"
#include "tools/Pbc.h"

#include <string>
#include <cmath>

//+PLUMEDOC MCOLVAR PLANES
/*
Calculate the components of the normal to the plane containing three atoms

\par Examples

*/
//+ENDPLUMEDOC

namespace PLMD {
namespace multicolvar {

class PlaneShortcut : public ActionShortcut {
private:
  void createVectorNormInput( const std::string& ilab, const std::string& olab, const std::string& vlab );
public:
  static void registerKeywords( Keywords& keys );
  PlaneShortcut(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(PlaneShortcut,"PLANES")

void PlaneShortcut::registerKeywords( Keywords& keys ) {
  ActionShortcut::registerKeywords( keys );
  MultiColvarShortcuts::shortcutKeywords( keys );
  keys.add("numbered","ATOMS","the sets of atoms that you would like to calculate the planes for");
  keys.add("numbered","LOCATION","the location at which the CV is assumed to be in space");
  keys.reset_style("LOCATION","atoms");
  keys.addFlag("VMEAN",false,"calculate the norm of the mean vector.");
  keys.addOutputComponent("_vmean","VMEAN","scalar","the norm of the mean vector");
  keys.addFlag("VSUM",false,"calculate the norm of the sum of all the vectors");
  keys.addOutputComponent("_vsum","VSUM","scalar","the norm of the mean vector");
  keys.needsAction("CENTER");
  keys.needsAction("GROUP");
  keys.needsAction("PLANE");
  keys.needsAction("MEAN");
  keys.needsAction("SUM");
  keys.needsAction("COMBINE");
  keys.needsAction("CUSTOM");
}

PlaneShortcut::PlaneShortcut(const ActionOptions&ao):
  Action(ao),
  ActionShortcut(ao) {
  bool vmean, vsum;
  parseFlag("VMEAN",vmean);
  parseFlag("VSUM",vsum);
  std::string dline;
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
  readInputLine( getShortcutLabel() + ": PLANE " + dline + " " + convertInputLineToString() );
  if( vmean ) {
    readInputLine( getShortcutLabel() + "_xs: MEAN ARG=" + getShortcutLabel() + ".x PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_ys: MEAN ARG=" + getShortcutLabel() + ".y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_zs: MEAN ARG=" + getShortcutLabel() + ".z PERIODIC=NO");
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vmean", "s" );
  }
  if( vsum ) {
    readInputLine( getShortcutLabel() + "_xz: SUM ARG=" + getShortcutLabel() + ".x PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_yz: SUM ARG=" + getShortcutLabel() + ".y PERIODIC=NO");
    readInputLine( getShortcutLabel() + "_zz: SUM ARG=" + getShortcutLabel() + ".z PERIODIC=NO");
    // Now calculate the total length of the vector
    createVectorNormInput( getShortcutLabel(), getShortcutLabel() + "_vsum", "z" );
  }
}

void PlaneShortcut::createVectorNormInput( const std::string& ilab, const std::string& olab, const std::string& vlab ) {
  readInputLine( olab + "2: COMBINE ARG=" + ilab + "_x" + vlab + "," + ilab + "_y" + vlab + "," + ilab + "_z" + vlab + " POWERS=2,2,2 PERIODIC=NO");
  readInputLine( olab + ": CUSTOM ARG=" + olab + "2 FUNC=sqrt(x) PERIODIC=NO");
}

}
}



