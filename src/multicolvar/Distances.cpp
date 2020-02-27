/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2020 The plumed team
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
#include "core/PlumedMain.h"
#include "core/Atoms.h"
#include "MultiColvarBase.h"

namespace PLMD {
namespace multicolvar {

class Distances : public ActionShortcut {
public:
  static void registerKeywords(Keywords& keys);
  explicit Distances(const ActionOptions&);
};

PLUMED_REGISTER_ACTION(Distances,"DISTANCES")

void Distances::registerKeywords(Keywords& keys) {
  ActionShortcut::registerKeywords( keys );
  keys.add("numbered","ATOMS","the pairs of atoms that you would like to calculate the angles for");
  keys.addFlag("NUMERICAL_DERIVATIVES", false, "calculate the derivatives for these quantities numerically");
  keys.addFlag("COMPONENTS",false,"calculate the x, y and z components of the distance separately and store them as label.x, label.y and label.z");
  keys.addFlag("SCALED_COMPONENTS",false,"calculate the a, b and c scaled components of the distance separately and store them as label.a, label.b and label.c");
  keys.reset_style("ATOMS","atoms"); MultiColvarBase::shortcutKeywords( keys );
  keys.add("atoms","ORIGIN","calculate the distance of all the atoms specified using the ATOMS keyword from this point");
}

Distances::Distances(const ActionOptions& ao):
Action(ao),
ActionShortcut(ao)
{
  // Create distances
  std::string dline = getShortcutLabel() + ": DISTANCE";
  bool numder; parseFlag("NUMERICAL_DERIVATIVES",numder); if( numder ) dline += " NUMERICAL_DERIVATIVES";
  bool comp; parseFlag("COMPONENTS",comp); if( comp ) dline += " COMPONENTS";
  bool scomp; parseFlag("SCALED_COMPONENTS",scomp); if( scomp ) dline += " SCALED_COMPONENTS";
  // Parse origin 
  std::string num, ostr; parse("ORIGIN",ostr);
  if( ostr.length()>0 ) {
      // Parse atoms
      std::vector<std::string> afstr, astr; parseVector("ATOMS",astr); Tools::interpretRanges(astr);
      for(unsigned i=0;i<astr.size();++i) {
          if( plumed.getAtoms().getAllGroups().count(astr[i]) ){
              const auto m=plumed.getAtoms().getAllGroups().find(astr[i]);
              for(unsigned j=0;j<m->second.size();++j) {
                  std::string num; Tools::convert( m->second[j].serial(), num ); afstr.push_back( num );
              }
          } else afstr.push_back(astr[i]);
      }
      for(unsigned i=0;i<afstr.size();++i) { Tools::convert( i+1, num ); dline += " ATOMS" + num + "=" + ostr + "," + afstr[i]; }
  } else {
      for(unsigned i=1;;++i) {
          std::string atstring; parseNumbered("ATOMS",i,atstring);
          if( atstring.length()==0 ) break;
          std::string num; Tools::convert( i, num );
          dline += " ATOMS" + num + "=" + atstring;
      }
  }
  readInputLine( dline );
  // Add shortcuts to label
  MultiColvarBase::expandFunctions( getShortcutLabel(), getShortcutLabel(), "", this );
}

}
}
